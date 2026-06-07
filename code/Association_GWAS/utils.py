#!/usr/bin/env python

import functools
import operator
import os
import re
import subprocess
import sys
import traceback
from io import StringIO
from itertools import groupby

import numpy as np
import pandas as pd
try:
    from natsort import natsorted
except ImportError:
    def natsorted(values):
        def key(value):
            parts = re.split(r"(\d+)", str(value))
            return [int(part) if part.isdigit() else part for part in parts]
        return sorted(values, key=key)


class NoTargetDataError(Exception):
    pass


def error_handler(func):
    @functools.wraps(func)
    def wrapper(num, *args, **kwargs):
        try:
            return func(num, *args, **kwargs)
        except NoTargetDataError:
            return None
        except Exception as e:
            e_info = sys.exc_info()
            e_type = e_info[0].__name__
            e_tracebk = ''.join(traceback.format_tb(e_info[2])[1:])
            print('Caught an exception for num={}:\n  {}: {}\nTraceback:\n{}'.format(num, e_type, e, e_tracebk))
        finally:
            sys.stdout.flush()

    return wrapper


def filter_vcf_line(line: bytes, bformat, col_inds, split_multi_GT):
    try:
        row = line.split(b'\t')
        data_fmts = row[8].split(b':')
        alt_alleles = row[4].split(b',')
        row[4] = alt_alleles[0]
        if len(data_fmts) > 1:
            data_ind = data_fmts.index(bformat)
            row[8] = bformat
            row = [row[x] if x <= 8 else row[x].split(b':')[data_ind] for x in col_inds]
        else:
            row = [row[x] for x in col_inds]
        if split_multi_GT and (len(alt_alleles) > 1):
            sample_str = b'\t'.join(row[4:])
            line_out = bytearray()
            for j in range(1, len(alt_alleles) + 1):
                str_j = sample_str
                for k in range(1, len(alt_alleles) + 1):
                    str_j = re.sub(str(k).encode(), b'.', str_j) if (k != j) else str_j
                str_j = re.sub(str(j).encode(), b'1', str_j)
                line_j = b'\t'.join([*row[0:3], alt_alleles[j - 1], str_j])
                line_j = line_j if line_j.endswith(b'\n') else line_j + b'\n'
                line_out += line_j
            return line_out
        line_out = b'\t'.join(row)
        line_out += b'' if line_out.endswith(b'\n') else b'\n'
        return line_out
    except Exception:
        return b''


def filter_weight_line(line: bytes, btarget: bytes, target_ind, col_inds):
    row = line.split(b'\t')
    if row[target_ind].startswith(btarget):
        line_out = b'\t'.join([row[x] for x in col_inds])
        line_out += b'' if line_out.endswith(b'\n') else b'\n'
        return line_out
    return b''


def filter_other_line(line: bytes, col_inds):
    row = line.split(b'\t')
    line_out = b'\t'.join([row[x] for x in col_inds])
    line_out += b'' if line_out.endswith(b'\n') else b'\n'
    return line_out


def get_snpIDs(df: pd.DataFrame, flip=False):
    chrms = df['CHROM'].astype('str').values
    pos = df['POS'].astype('str').values
    ref = df['REF'].values
    alt = df['ALT'].values
    if flip:
        return [':'.join(i) for i in zip(chrms, pos, alt, ref)]
    return [':'.join(i) for i in zip(chrms, pos, ref, alt)]


def reformat_sample_vals(df: pd.DataFrame, data_format, sampleID):
    df = df.reset_index(drop=True).copy()
    vals = df[sampleID].values
    if data_format == 'GT':
        vals[(vals == '0|0') | (vals == '0/0')] = 0
        vals[(vals == '1|0') | (vals == '1/0') | (vals == '0|1') | (vals == '0/1')] = 1
        vals[(vals == '1|1') | (vals == '1/1')] = 2
        vals[(vals == '.|.') | (vals == './.')] = np.nan
    elif data_format == 'DS':
        vals[(vals == '.')] = np.nan
    vals = vals.astype(np.float32)
    df = pd.concat([df.drop(columns=sampleID), pd.DataFrame(vals, columns=sampleID)], axis=1)
    return df


def read_tabix(start, end, sampleID, chrm, path, file_cols, col_inds, cols, dtype,
               genofile_type=None, data_format=None, target_ind=5, target=None,
               weight_threshold=0, raise_error=True, **kwargs):
    command_str = ' '.join(['tabix', path, chrm + ':' + start + '-' + end])
    proc = subprocess.Popen([command_str], shell=True, stdout=subprocess.PIPE)
    proc_out = bytearray()

    if genofile_type == 'vcf':
        bformat = str.encode(data_format)
        filter_line = functools.partial(filter_vcf_line, bformat=bformat, col_inds=col_inds,
                                        split_multi_GT=data_format == 'GT')
    elif genofile_type == 'weight':
        btarget = str.encode(target)
        filter_line = functools.partial(filter_weight_line, btarget=btarget, target_ind=target_ind, col_inds=col_inds)
    elif genofile_type == 'bgw_weight':
        btarget = b'0'
        filter_line = functools.partial(filter_weight_line, btarget=btarget, target_ind=target_ind, col_inds=col_inds)
    else:
        filter_line = functools.partial(filter_other_line, col_inds=col_inds)

    while proc.poll() is None:
        line = proc.stdout.readline()
        if len(line) == 0:
            break
        proc_out += filter_line(line)
    for line in proc.stdout:
        proc_out += filter_line(line)

    if not proc_out and raise_error:
        print('No tabix data for target.')
        raise NoTargetDataError

    df = pd.read_csv(
        StringIO(proc_out.decode('utf-8')),
        sep='\t',
        low_memory=False,
        header=None,
        names=cols,
        dtype=dtype)

    if len(sampleID):
        df = df[df[sampleID].count(axis=1) != 0].reset_index(drop=True)

    df = optimize_cols(df)

    if (genofile_type != 'weight') or ('snpID' not in cols):
        df['snpID'] = get_snpIDs(df)

    df = df.drop_duplicates(['snpID'], keep='first').reset_index(drop=True)

    if genofile_type == 'weight':
        if 'ES' not in cols:
            if ('b' in cols) and ('beta' in cols):
                df['ES'] = df['b'] + df['beta']
            if ('PCP' in cols) and ('beta' in cols):
                df['ES'] = np.prod(df['PCP'], df['beta'])

        if weight_threshold:
            df = df[operator.gt(np.abs(df['ES']), weight_threshold)].reset_index(drop=True)
            if df.empty and raise_error:
                print('No test SNPs above weight threshold for target: ' + target + '.')
                raise NoTargetDataError

    if data_format == 'GT':
        valid_GT = ['.|.', '0|0', '0|1', '1|0', '1|1',
                    './.', '0/0', '0/1', '1/0', '1/1']
        df = df[np.all(df[sampleID].isin(valid_GT), axis=1)].reset_index(drop=True)

    if data_format in ('GT', 'DS'):
        df = reformat_sample_vals(df, data_format, sampleID)

    if df.empty and raise_error:
        print('No valid tabix data for target.')
        raise NoTargetDataError

    return df


def format_elapsed_time(time_secs):
    val = abs(int(time_secs))
    day = val // (3600 * 24)
    hour = val % (3600 * 24) // 3600
    mins = val % 3600 // 60
    secs = val % 60
    res = '%02d:%02d:%02d:%02d' % (day, hour, mins, secs)
    if int(time_secs) < 0:
        res = "-%s" % res
    return res


def get_header(path, out='tuple', zipped=False, rename={}):
    compress_type = 'gzip' if zipped else None
    rename = {**{'#CHROM': 'CHROM'}, **rename}
    header = pd.read_csv(
        path,
        sep='\t',
        header=0,
        compression=compress_type,
        low_memory=False,
        nrows=0).rename(columns=rename)
    if out == 'tuple':
        return tuple(header)
    if out == 'list':
        return list(header)
    return header


def get_cols_dtype(file_cols, cols, sampleid=None, genofile_type=None, add_cols=[], drop_cols=[], get_id=False,
                   ind_namekey=False, **kwargs):
    if sampleid is not None:
        cols = cols + sampleid.tolist()
        sampleid_dtype = object if genofile_type == 'vcf' else np.float64
        sampleid_dict = {x: sampleid_dtype for x in sampleid}
    else:
        sampleid_dict = {}

    dtype_dict = {
        'ALT': object,
        'b': np.float64,
        'beta': np.float64,
        'BETA': np.float64,
        'CHROM': object,
        'COV': object,
        'ES': np.float64,
        'FILTER': object,
        'FORMAT': object,
        'INFO': object,
        'GeneEnd': np.int64,
        'GeneName': object,
        'GeneStart': np.int64,
        'ID': object,
        'MAF': np.float64,
        'PCP': np.float64,
        'POS': np.int64,
        'QUAL': object,
        'REF': object,
        'SE': np.float64,
        'snpID': object,
        'TargetID': object,
        'Trans': np.int64,
        'Zscore': np.float64,
        **sampleid_dict}

    cols = cols + add_cols
    if get_id:
        if 'snpID' in file_cols:
            cols.append('snpID')
        elif 'ID' in file_cols:
            cols.append('ID')

    cols = [x for x in cols if (x not in drop_cols)]
    col_inds = tuple(sorted([file_cols.index(x) for x in cols]))

    if ind_namekey:
        out_dtype_dict = {x: dtype_dict[x] for x in cols}
    else:
        ind_dtype_dict = {file_cols.index(x): dtype_dict[x] for x in cols}
        out_dtype_dict = ind_dtype_dict

    return {
        'file_cols': file_cols,
        'cols': [file_cols[i] for i in col_inds],
        'col_inds': col_inds,
        'dtype': out_dtype_dict}


def get_ld_cols(path):
    file_cols = tuple(pd.read_csv(
        path,
        sep='\t',
        header=0,
        compression='gzip',
        low_memory=False,
        nrows=0).rename(
        columns={'#snpID': 'snpID', '#ID': 'snpID', 'ID': 'snpID', '#0': 'row', '#row': 'row', '0': 'row'}))
    cols = ['row', 'snpID', 'COV']
    file_cols_ind = tuple([file_cols.index(x) for x in cols])
    return file_cols, file_cols_ind


def get_ld_regions_list(snp_ids):
    chrm = snp_ids[0].split(':')[0] + ':'
    pos_vals = [int(snp.split(':')[1]) for snp in snp_ids]
    pos_vals = sorted(list(set(pos_vals)))
    for x, y in groupby(enumerate(pos_vals), lambda p: p[1] - p[0]):
        y = list(y)
        yield chrm + str(y[0][1]) + '-' + str(y[-1][1])


def call_tabix_regions(path, regs_str, filter_line=lambda x: x):
    commond = ['tabix ' + path + ' ' + regs_str]
    proc = subprocess.Popen(commond, shell=True, stdout=subprocess.PIPE)
    proc_out = bytearray()
    while proc.poll() is None:
        line = proc.stdout.readline()
        if len(line) == 0:
            break
        proc_out += filter_line(line)
    for line in proc.stdout:
        proc_out += filter_line(line)
    return proc_out


def get_ld_regions_data(regs_str, path, snp_ids, ld_cols, ld_cols_ind):
    proc_out = call_tabix_regions(path, regs_str)
    regs_data = pd.read_csv(
        StringIO(proc_out.decode('utf-8')),
        sep='\t',
        low_memory=False,
        header=None,
        names=ld_cols,
        usecols=ld_cols_ind,
        dtype={'snpID': object, 'row': np.int32, 'COV': object}
    ).drop_duplicates(['snpID'], keep='first')
    regs_data = regs_data[regs_data.snpID.isin(snp_ids)]
    return regs_data


def get_ld_data(path, snp_ids):
    ld_cols, ld_cols_ind = get_ld_cols(path)
    regs_lst = list(get_ld_regions_list(snp_ids))
    N = len(regs_lst)
    regs_args = [path, snp_ids, ld_cols, ld_cols_ind]
    try:
        regs_str = ' '.join(regs_lst)
        cov_data = get_ld_regions_data(regs_str, *regs_args)
    except OSError:
        n = 2500
        while n:
            try:
                regs_str_lst = [' '.join(regs_lst[i:i + n]) for i in range(0, N, n)]
                cov_data = pd.concat([get_ld_regions_data(regs_str, *regs_args) for regs_str in regs_str_lst])
            except OSError:
                n -= 500
                pass
            else:
                n = 0
    return cov_data.set_index('row')


def get_ld_matrix(MCOV, return_diag=False):
    MCOV = MCOV.copy()
    MCOV['COV'] = MCOV['COV'].apply(lambda x: np.fromstring(x, dtype=np.float32, sep=','))
    inds = MCOV.index
    n_inds = inds.size
    V_upper = np.zeros((n_inds, n_inds))
    for i in range(n_inds):
        cov_i = MCOV.COV.loc[inds[i]]
        N = cov_i.size
        for j in range(i, n_inds):
            if inds[j] - inds[i] < N:
                V_upper[i, j] = cov_i[inds[j] - inds[i]]
            else:
                V_upper[i, j] = 0
    snp_Var = V_upper.diagonal()
    V = V_upper + V_upper.T - np.diag(snp_Var)
    snp_sd = np.sqrt(snp_Var)
    if return_diag:
        return snp_sd, V, snp_Var
    return snp_sd, V


complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def get_complement_snpID(snpID):
    if not isinstance(snpID, str):
        return None
    try:
        parts = snpID.split(':')
        if len(parts) != 4:
            return None
        chrom, pos, ref, alt = parts
        if not ref or not alt:
            return None
        ref_comp_chars = []
        for base in ref.upper():
            comp = complement_map.get(base)
            if comp is None:
                return None
            ref_comp_chars.append(comp)
        ref_comp = "".join(ref_comp_chars)
        alt_comp_chars = []
        for base in alt.upper():
            comp = complement_map.get(base)
            if comp is None:
                return None
            alt_comp_chars.append(comp)
        alt_comp = "".join(alt_comp_chars)
        return f"{chrom}:{pos}:{ref_comp}:{alt_comp}"
    except Exception:
        return None


def get_complement_snpID_2(snpID):
    if not isinstance(snpID, str):
        return None
    try:
        parts = snpID.rsplit(':', 3)
        if len(parts) != 4:
            return None
        chrom, pos, ref, alt = parts
        if not ref or not alt:
            return None
        ref_comp_chars = []
        for base in ref.upper():
            comp = complement_map.get(base)
            if comp is None:
                ref_comp_chars.append(base)
            else:
                ref_comp_chars.append(comp)
        ref_comp = "".join(ref_comp_chars)
        alt_comp_chars = []
        for base in alt.upper():
            comp = complement_map.get(base)
            if comp is None:
                alt_comp_chars.append(base)
            else:
                alt_comp_chars.append(comp)
        alt_comp = "".join(alt_comp_chars)
        return f"{chrom}:{pos}:{ref_comp}:{alt_comp}"
    except Exception:
        return None


def flip_snpIDs(snpIDs):
    return np.array([':'.join([y[0], y[1], y[3], y[2]]) for y in [x.split(':') for x in snpIDs]])


def optimize_cols(df: pd.DataFrame):
    if 'CHROM' in df.columns:
        try:
            df['CHROM'] = df['CHROM'].astype(str).astype(np.int8)
        except ValueError:
            pass
    ints = df.select_dtypes(include=['int64']).columns.tolist()
    df[ints] = df[ints].apply(pd.to_numeric, downcast='integer')
    floats = df.select_dtypes(include=['float64']).columns.tolist()
    df[floats] = df[floats].apply(pd.to_numeric, downcast='float')
    return df


def merge_regions(regions):
    if not regions:
        return []
    df = pd.DataFrame(regions, columns=['CHROM', 'Start', 'End'])
    df['Start'] = pd.to_numeric(df['Start'])
    df['End'] = pd.to_numeric(df['End'])
    df['CHROM'] = df['CHROM'].astype(str)

    df['CHROM_SortKey'] = df['CHROM']
    df = df.sort_values(by=['CHROM_SortKey', 'Start'], key=lambda x: x if x.name == 'Start' else natsorted(x))
    df = df.drop(columns=['CHROM_SortKey'])

    merged = []
    if df.empty:
        return merged

    for chrom, group in df.groupby('CHROM', sort=False):
        if group.empty:
            continue
        group_list = group[['Start', 'End']].values.tolist()
        if not group_list:
            continue
        current_start, current_end = group_list[0]
        for i in range(1, len(group_list)):
            next_start, next_end = group_list[i]
            if next_start <= current_end + 1:
                current_end = max(current_end, next_end)
            else:
                merged.append((str(chrom), current_start, current_end))
                current_start, current_end = next_start, next_end
        merged.append((str(chrom), current_start, current_end))
    return merged
