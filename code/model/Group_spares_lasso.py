import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error
from group_lasso import GroupLasso
import statsmodels.api as sm
from sklearn.linear_model import ElasticNetCV, LassoCV, RidgeCV
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.metrics import make_scorer, mean_squared_error
from scipy.stats import uniform, pearsonr
import os
import tempfile
import tempfile
import sys
import subprocess

# Add parent directory to path for importing GRNutils
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import GRNutils as tg

def mse_scorer(estimator, X, y):
    y_pred = estimator.predict(X)
    return mean_squared_error(y, y_pred)


# Compatibility wrapper used by downstream wrappers
def sparse_group_lasso_Qualified(groups, train_geno_exp, test_geno_exp):
    return compare_lasso_enet_cv(groups, train_geno_exp, test_geno_exp)


def run_dpr(Bimbam_df, Pheno_df, output_dir, target_id, dpr_path, method, es_type):
    if not dpr_path or (not os.path.exists(dpr_path)):
        print(f"： DPR ：{dpr_path}")
        return None

    bimbam_file = f"{output_dir}{target_id}.bimbam"
    pheno_file = f"{output_dir}{target_id}.pheno"
    out_args = {'header': False, 'index': None, 'sep': '\t', 'mode': 'w', 'float_format': '%f'}
    Bimbam_df.to_csv(bimbam_file, **out_args)
    Pheno_df.to_csv(pheno_file, **out_args)
    output_prefix = f"{output_dir}output/{target_id}.param.txt"
    cmd = [dpr_path, '-g', f"{target_id}.bimbam", '-p', f"{target_id}.pheno", '-dpr', method,
           '-notsnp',
           '-o', target_id]
    print(cmd)
    try:
        subprocess.check_call(cmd, cwd=output_dir)
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"：DPR ， DPR 。: {e}")
        return None
    finally:
        os.remove(bimbam_file)
        os.remove(pheno_file)

    print("esx")  # 

    if os.path.exists(output_prefix):
        dpr_result = pd.read_csv(output_prefix,
                                 sep='\t',
                                 header=0,
                                 names=['CHROM', 'snpID', 'POS', 'n_miss', 'b', 'beta', 'gamma'],
                                 usecols=['snpID', 'b', 'beta'],
                                 dtype={'snpID': object, 'b': np.float64, 'beta': np.float64})
        os.remove(output_prefix)
        os.remove(f"{output_dir}output/{target_id}.log.txt")
        DPR_Out = tg.optimize_cols(dpr_result)

        # GET EFFECT SIZE
        DPR_Out['ES'] = DPR_Out['beta']

        return DPR_Out

    else:
        print(f"：DPR  {output_prefix} ")
        return None


def compare_lasso_enet_cv(groups, train, test=None, k=5, random_state=42, Geno_meta=None, dpr_path=None, dpr_method='1',
                          es_type='fixed', tmp_DPR=None, used_model=None):

    trainX = train.iloc[:, :-1]
    trainY = train.iloc[:, -1].values
    if test is not None:
        testX = test.iloc[:, :-1]
        testY = test.iloc[:, -1].values
    else:
        testX = trainX
        testY = trainY

    if groups is not None:
        gl = GroupLasso(
            groups=groups, l1_reg=0, group_reg=0.1, scale_reg="inverse_group_size",
            supress_warning=True, frobenius_lipschitz=True, fit_intercept=False
        )
        gl.fit(trainX.values, trainY)
        nonzero_groups = [i for i, group in enumerate(gl.groups_) if np.any(gl.coef_[group] != 0)]
        selected_features = [trainX.columns[gl.groups_[group_index]] for group_index in nonzero_groups]
        selected_features = [item for sublist in selected_features for item in sublist]
    else:
        selected_features = trainX.columns.tolist()
    if not selected_features:
        print("：Group Lasso 。。")
        selected_features = trainX.columns.tolist()

    # LassoCV
    lasso_cv_model = LassoCV(
        cv=k, n_jobs=max(1, int(0.7 * os.cpu_count())), random_state=random_state,
        fit_intercept=True, alphas=np.array([0.0001, 0.001, 0.01, 0.05, 0.1, 0.15, 0.3, 0.5])
    )
    print(" LassoCV ...")
    lasso_cv_model.fit(trainX[selected_features].values, trainY)
    lasso_predY = lasso_cv_model.predict(testX[selected_features].values)
    lasso_lm = sm.OLS(testY, sm.add_constant(lasso_predY)).fit()
    lasso_Rsquared = lasso_lm.rsquared
    print("LassoCV 。")

    # ElasticNetCV
    enet_cv_model = ElasticNetCV(
        max_iter=100, l1_ratio=[0.05, 0.1, 0.2, 0.5], fit_intercept=True,
        alphas=np.array([0.0001, 0.001, 0.01, 0.05, 0.15]), selection='random',
        cv=k, n_jobs=max(1, int(0.7 * os.cpu_count())), random_state=random_state
    )
    print(" ElasticNetCV ...")
    enet_cv_model.fit(trainX[selected_features].values, trainY)
    enet_predY = enet_cv_model.predict(testX[selected_features].values)
    enet_lm = sm.OLS(testY, sm.add_constant(enet_predY)).fit()
    enet_Rsquared = enet_lm.rsquared
    print("ElasticNetCV 。")

    # DPR
    Rsquared_dpr = -np.inf
    predY_dpr = None
    lm_dpr = None
    weights_dpr = None
    if dpr_path is not None and Geno_meta is not None:
        print(" DPR ...")
        selected_snps = selected_features
        meta_info = Geno_meta[Geno_meta['snpID'].isin(selected_snps)][['snpID', 'REF', 'ALT']].copy()

        geno_data = trainX[selected_snps].T
        geno_data.index.name = 'snpID'

        Bimbam_df = meta_info.merge(geno_data, left_on='snpID', right_index=True)

        sample_columns = [col for col in Bimbam_df.columns if col not in ['snpID', 'REF', 'ALT']]
        Bimbam_df = Bimbam_df[['snpID', 'REF', 'ALT'] + sample_columns]

        Pheno_df = pd.DataFrame(trainY, columns=['expression'])

        dpr_output_dir = tmp_DPR if tmp_DPR else tempfile.mkdtemp(prefix='grntwas_dpr_')
        os.makedirs(dpr_output_dir, exist_ok=True)
        os.makedirs(f"{dpr_output_dir}output", exist_ok=True)
        target_id = train.columns[-1]
        dpr_result = run_dpr(Bimbam_df, Pheno_df, dpr_output_dir, target_id, dpr_path, dpr_method, es_type)
        if dpr_result is not None:
            weights_dpr = dpr_result.set_index('snpID')['ES']
            common_snps = [snp for snp in weights_dpr.index if snp in testX.columns]
            if common_snps:
                testX_selected = testX[common_snps]
                predY_dpr = np.dot(testX_selected.values, weights_dpr[common_snps].values)
                lm_dpr = sm.OLS(testY, sm.add_constant(predY_dpr)).fit()
                Rsquared_dpr = lm_dpr.rsquared
            else:
                print("：DPR  SNP ")
        else:
            print("：DPR ")
        print("DPR 。")
    models = {
        "LassoCV": {
            "Rsquared": lasso_Rsquared,
            "predY": lasso_predY,
            "model": lasso_cv_model,
            "lm": lasso_lm,
            "Alpha": None,  # LassoCV  l1_ratio
            "Lambda": lasso_cv_model.alpha_,
            "cvm": np.min(lasso_cv_model.mse_path_),
            "beta": lasso_cv_model.coef_
        },
        "ElasticNetCV": {
            "Rsquared": enet_Rsquared,
            "predY": enet_predY,
            "model": enet_cv_model,
            "lm": enet_lm,
            "Alpha": enet_cv_model.l1_ratio_,
            "Lambda": enet_cv_model.alpha_,
            "cvm": np.min(enet_cv_model.mse_path_),
            "beta": enet_cv_model.coef_
        },
        "DPR": {
            "Rsquared": Rsquared_dpr,
            "predY": predY_dpr,
            "model": None,  # DPR 
            "lm": lm_dpr,
            "Alpha": None,
            "Lambda": None,
            "cvm": None,
            "beta": weights_dpr  # DPR 
        }
    }
    if used_model:
        best_model_name = used_model
    else:
        best_model_name = max(models, key=lambda x: models[x]["Rsquared"])

    best_info = models[best_model_name]
    print(f": {best_model_name}, R² = {best_info['Rsquared']}")

    if test is not None:
        return best_info["Rsquared"], best_model_name
    else:
        beta_model = best_info["beta"]
        Rsquared = best_info["Rsquared"]
        Pvalue = best_info["lm"].f_pvalue if best_info["lm"] is not None else None
        Alpha = best_info["Alpha"]
        Lambda = best_info["Lambda"]
        cvm = best_info["cvm"]

        beta = np.zeros(trainX.shape[1])
        if isinstance(beta_model, pd.Series):  # DPR 
            for snp in beta_model.index:
                if snp in selected_features:
                    feature_index = trainX.columns.get_loc(snp)
                    beta[feature_index] = beta_model[snp]
        else:  # LassoCV  ElasticNetCV 
            for i, feature_name in enumerate(selected_features):

                feature_index = trainX.columns.get_loc(feature_name)
                beta[feature_index] = beta_model[i]

        return beta, Rsquared, Pvalue, Alpha, Lambda, cvm


def compare_lasso_enet_cv_revise(groups, train, test=None, k=5, random_state=42, Geno_meta=None, dpr_path=None, dpr_method='1',
                          es_type='fixed', tmp_DPR=None, used_model='dpr'):
    trainX = train.iloc[:, :-1]
    trainY = train.iloc[:, -1].values
    if test is not None:
        testX = test.iloc[:, :-1]
        testY = test.iloc[:, -1].values
    else:
        testX = trainX
        testY = trainY

    if groups is not None:
        gl = GroupLasso(
            groups=groups, l1_reg=0, group_reg=0.1, scale_reg="inverse_group_size",
            supress_warning=True, frobenius_lipschitz=True, fit_intercept=False
        )
        gl.fit(trainX.values, trainY)
        nonzero_groups = [i for i, group in enumerate(gl.groups_) if np.any(gl.coef_[group] != 0)]
        selected_features = [trainX.columns[gl.groups_[group_index]] for group_index in nonzero_groups]
        selected_features = [item for sublist in selected_features for item in sublist]
    else:
        selected_features = trainX.columns.tolist()
    if not selected_features:
        print("GroupLasso selected no features; using all SNPs.")
        selected_features = trainX.columns.tolist()

    requested_model = (used_model or 'dpr').lower()
    valid_models = {'dpr', 'elasticnet', 'lasso', 'ridge', 'auto'}
    if requested_model not in valid_models:
        raise ValueError(f"Unsupported model '{used_model}'. Choose from {sorted(valid_models)}.")

    model_names = ['dpr', 'elasticnet', 'lasso', 'ridge'] if requested_model == 'auto' else [requested_model]
    models = {}

    def add_prediction_model(name, pred_y, model, alpha=None, lam=None, cvm=None, beta=None):
        lm = sm.OLS(testY, sm.add_constant(pred_y)).fit()
        models[name] = {
            "Rsquared": lm.rsquared,
            "predY": pred_y,
            "model": model,
            "lm": lm,
            "Alpha": alpha,
            "Lambda": lam,
            "cvm": cvm,
            "beta": beta
        }

    if 'dpr' in model_names:
        if dpr_path is None or not os.path.exists(dpr_path):
            if requested_model == 'dpr':
                raise FileNotFoundError(
                    "DPR is the default GRNTWAS training model, but the DPR executable was not found. "
                    "Provide --dpr-path /path/to/DPR, or explicitly use --model elasticnet, --model lasso, --model ridge, or --model auto."
                )
            print("DPR executable not found; skipping DPR for auto model selection.")
        elif Geno_meta is None:
            if requested_model == 'dpr':
                raise ValueError("DPR requires Geno_meta SNP metadata, but Geno_meta was not provided.")
            print("Geno_meta not provided; skipping DPR for auto model selection.")
        else:
            print(" DPR ...")
            selected_snps = selected_features
            meta_info = Geno_meta[Geno_meta['snpID'].isin(selected_snps)][['snpID', 'REF', 'ALT']].copy()
            geno_data = trainX[selected_snps].T
            geno_data.index.name = 'snpID'
            Bimbam_df = meta_info.merge(geno_data, left_on='snpID', right_index=True)
            sample_columns = [col for col in Bimbam_df.columns if col not in ['snpID', 'REF', 'ALT']]
            Bimbam_df = Bimbam_df[['snpID', 'REF', 'ALT'] + sample_columns]
            Pheno_df = pd.DataFrame(trainY, columns=['expression'])
            dpr_output_dir = tmp_DPR if tmp_DPR else tempfile.mkdtemp(prefix='grntwas_dpr_')
            os.makedirs(dpr_output_dir, exist_ok=True)
            os.makedirs(f"{dpr_output_dir}output", exist_ok=True)
            target_id = train.columns[-1]
            dpr_result = run_dpr(Bimbam_df, Pheno_df, dpr_output_dir, target_id, dpr_path, dpr_method, es_type)
            if dpr_result is not None:
                weights_dpr = dpr_result.set_index('snpID')['ES']
                common_snps = [snp for snp in weights_dpr.index if snp in testX.columns]
                if common_snps:
                    testX_selected = testX[common_snps]
                    predY_dpr = np.dot(testX_selected.values, weights_dpr[common_snps].values)
                    lm_dpr = sm.OLS(testY, sm.add_constant(predY_dpr)).fit()
                    models["DPR"] = {
                        "Rsquared": lm_dpr.rsquared,
                        "predY": predY_dpr,
                        "model": None,
                        "lm": lm_dpr,
                        "Alpha": None,
                        "Lambda": None,
                        "cvm": None,
                        "beta": weights_dpr
                    }
                else:
                    print("DPR returned no SNPs present in the test matrix.")
            elif requested_model == 'dpr':
                raise RuntimeError("DPR ran but did not produce a usable SNP-weight result.")
            print("DPR finished.")

    if 'elasticnet' in model_names:
        print(" ElasticNetCV ...")
        enet_cv_model = ElasticNetCV(cv=k, n_jobs=max(1, int(0.7 * os.cpu_count())), random_state=random_state)
        enet_cv_model.fit(trainX[selected_features].values, trainY)
        enet_predY = enet_cv_model.predict(testX[selected_features].values)
        add_prediction_model(
            "ElasticNetCV",
            enet_predY,
            enet_cv_model,
            alpha=enet_cv_model.l1_ratio_,
            lam=enet_cv_model.alpha_,
            cvm=np.min(enet_cv_model.mse_path_),
            beta=enet_cv_model.coef_
        )
        print("ElasticNetCV finished.")

    if 'lasso' in model_names:
        print(" LassoCV ...")
        lasso_cv_model = LassoCV(cv=k, n_jobs=max(1, int(0.7 * os.cpu_count())), random_state=random_state)
        lasso_cv_model.fit(trainX[selected_features].values, trainY)
        lasso_predY = lasso_cv_model.predict(testX[selected_features].values)
        add_prediction_model(
            "LassoCV",
            lasso_predY,
            lasso_cv_model,
            alpha=None,
            lam=lasso_cv_model.alpha_,
            cvm=np.min(lasso_cv_model.mse_path_),
            beta=lasso_cv_model.coef_
        )
        print("LassoCV finished.")

    if 'ridge' in model_names:
        print(" RidgeCV ...")
        ridge_cv_model = RidgeCV()
        ridge_cv_model.fit(trainX[selected_features].values, trainY)
        ridge_predY = ridge_cv_model.predict(testX[selected_features].values)
        add_prediction_model(
            "RidgeCV",
            ridge_predY,
            ridge_cv_model,
            alpha=None,
            lam=getattr(ridge_cv_model, 'alpha_', None),
            cvm=None,
            beta=ridge_cv_model.coef_
        )
        print("RidgeCV finished.")

    if not models:
        raise RuntimeError(f"No usable model was trained for requested model '{requested_model}'.")

    best_model_name = max(models, key=lambda x: models[x]["Rsquared"])
    best_info = models[best_model_name]
    print(f": {best_model_name}, R² = {best_info['Rsquared']}")

    beta_model = best_info["beta"]
    if beta_model is None:
        raise RuntimeError(f"Selected model {best_model_name} did not produce SNP weights.")

    if test is not None:
        return best_info["Rsquared"], best_model_name

    Rsquared = best_info["Rsquared"]
    Pvalue = best_info["lm"].f_pvalue if best_info["lm"] is not None else None
    Alpha = best_info["Alpha"]
    Lambda = best_info["Lambda"]
    cvm = best_info["cvm"]

    beta = np.zeros(trainX.shape[1])
    if isinstance(beta_model, pd.Series):
        for snp in beta_model.index:
            if snp in trainX.columns:
                feature_index = trainX.columns.get_loc(snp)
                beta[feature_index] = beta_model[snp]
    else:
        for i, feature_name in enumerate(selected_features):
            feature_index = trainX.columns.get_loc(feature_name)
            beta[feature_index] = beta_model[i]

    return beta, Rsquared, Pvalue, Alpha, Lambda, cvm, best_model_name
