import os
import pandas as pd
import numpy as np
import argparse
import warnings
import random
import re
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, NuSVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, accuracy_score, roc_auc_score, precision_recall_curve, auc, confusion_matrix
from xgboost import XGBClassifier

# Suppress warnings
warnings.filterwarnings('ignore')

# Configuration
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "Dataset")
RESULTS_PP_DIR = os.path.join(BASE_DIR, "Results_PP")
FS_RESULTS_PATH = os.path.join(BASE_DIR, "filtered_FS_results.csv")
os.makedirs(RESULTS_PP_DIR, exist_ok=True)

# Full Configuration Lists
ALL_VIRUSES = ["H1N1", "H3N2", "HRV", "RSV"]
ALL_LEVELS = {
    "Probe": "probe_expression", 
    "ssGSEA": "gene_ssgsea", 
    "Combined": "combined"
}

# CSV Feature Mapping (Script Name -> CSV Name)
CSV_FEATURE_MAP = {
    "Probe": "probe_expression",
    "ssGSEA": "gene_ssgsea",
    "Combined": "probe_and_ssgsea"
}

ALL_CLASSIFIERS = ["LR", "XGB", "SVM", "KNN", "NuSVC", "GNB", "RF"]
TIME_POINTS = [0, 24, 48, 72, 96, 120]

# Requested Experiment Order
EXPERIMENT_ORDER = [
    "RSV_DEE1", 
    "H3N2_DEE2", 
    "H1N1_DEE3", 
    "H1N1_DEE4", 
    "H3N2_DEE5", 
    "HRV_DUKE", 
    "HRV_UVA"
]

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)

def get_model(classifier_name, seed):
    if classifier_name == "LR":
        return LogisticRegression(max_iter=1000000,solver="liblinear")
    elif classifier_name == "XGB":
        return XGBClassifier(random_state=seed)
    elif classifier_name == "SVM":
        return SVC(probability=True, random_state=seed)
    elif classifier_name == "NuSVC":
        return NuSVC(probability=True, random_state=seed)
    elif classifier_name == "KNN":
        return KNeighborsClassifier()
    elif classifier_name == "GNB":
        return GaussianNB()
    elif classifier_name == "RF":
        return RandomForestClassifier(random_state=seed)
    else:
        raise ValueError(f"Unknown classifier: {classifier_name}")

def train_predict(model, X_train, y_train, X_test):
    model.fit(X_train, y_train)
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X_test)[:, 1]
    else:
        return model.predict(X_test)

def get_sort_key(subject_id):
    """
    Returns a tuple (experiment_index, subject_number) for sorting.
    """
    s_id = str(subject_id)
    exp_index = float('inf')
    for i, exp in enumerate(EXPERIMENT_ORDER):
        if s_id.startswith(exp):
            exp_index = i
            break
    nums = [int(c) for c in re.split(r'(\d+)', s_id) if c.isdigit()]
    subj_num = nums[-1] if nums else 0
    return (exp_index, subj_num)

def load_fs_data():
    if os.path.exists(FS_RESULTS_PATH):
        return pd.read_csv(FS_RESULTS_PATH)
    else:
        print(f"Warning: Feature Selection file not found at {FS_RESULTS_PATH}")
        return None

def get_selected_features(fs_df, fs_method, virus, classifier, target, tp, feature_level):
    """
    Filters the FS DataFrame to find the selected features for the specific combination.
    """
    if fs_df is None:
        return None
        
    # Map script names to CSV names
    csv_feature = CSV_FEATURE_MAP.get(feature_level)
    csv_classifier = classifier + "*" # CSV uses "LR*" etc.
    csv_sc = int(target[-1]) # "SC1" -> 1
    
    # Filter
    subset = fs_df[
        (fs_df["Experiment"] == virus) &
        (fs_df["Classifier"] == csv_classifier) &
        (fs_df["SC"] == csv_sc) &
        (fs_df["TP"] == tp) &
        (fs_df["Feature"] == csv_feature) & 
        (fs_df["FS Method"] == fs_method)
    ]
    
    if subset.empty:
        return None
        
    feat_str = subset.iloc[0]["SelectedFeatures"]
    if pd.isna(feat_str) or feat_str == "":
        return []
        
    return str(feat_str).split(";")

def run_analysis(classifiers, targets, viruses, features, time_points, seed, fs_method=None):
    print(f"Starting Analysis with Seed={seed}")
    print(f"  Classifiers: {classifiers}")
    print(f"  Targets: {targets}")
    print(f"  Viruses: {viruses}")
    print(f"  Features: {list(features.keys())}")
    print(f"  TimePoints: {time_points}")
    print(f"  FS Method: {fs_method if fs_method else 'None'}")
    
    set_seed(seed)
    
    # Load FS Data if needed
    fs_df = None
    if fs_method:
        fs_df = load_fs_data()
        if fs_df is None:
            print("  Error: FS Method requested but FS file missing/empty. Aborting.")
            return

    summary_metrics = []
    
    # Iterate Filters
    for classifier_name in classifiers:
        for target in targets:
            for level_name, level_dir in features.items():
                print(f"Processing: {classifier_name} | {target} | {level_name}")
                
                for t_time in time_points:
                    print(f"  T{t_time}")
                    aggregated_predictions = []
                    
                    for virus in viruses:
                        
                        # FS Lookup for THIS virus
                        current_selected_feats = None
                        if fs_method:
                            current_selected_feats = get_selected_features(fs_df, fs_method, virus, classifier_name, target, t_time, level_name)
                            if current_selected_feats is None:
                                print(f"    Warning: No FS entry found for {virus} {fs_method}. Skipping.")
                                continue
                            if len(current_selected_feats) == 0:
                                print(f"    Warning: 0 features selected for {virus} {fs_method}. Skipping.")
                                continue

                        # Construct file paths
                        if level_name == "Probe": file_type = "Probe"
                        elif level_name == "ssGSEA": file_type = "ssGSEA"
                        else: file_type = "Combined"
                        
                        train_file = os.path.join(DATA_DIR, virus, level_dir, f"{virus}_{file_type}_train_T{t_time}.csv")
                        test_file = os.path.join(DATA_DIR, virus, level_dir, f"{virus}_{file_type}_test_T{t_time}.csv")
                        
                        if not os.path.exists(train_file) or not os.path.exists(test_file):
                            continue
                        
                        # Load Data
                        try:
                            train_df = pd.read_csv(train_file)
                            test_df = pd.read_csv(test_file)
                        except Exception as e:
                            print(f"    Error reading file for {virus}: {e}")
                            continue
                        
                        if target not in train_df.columns:
                            continue
                        
                        # Prepare Data
                        drop_cols = ["SubjectID", "SC1", "SC2"]
                        feature_cols = [c for c in train_df.columns if c not in drop_cols]
                        
                        # Apply Feature Selection if Active
                        if fs_method:
                            valid_feats = [f for f in current_selected_feats if f in feature_cols]
                            if len(valid_feats) == 0:
                                print(f"    Error: None of the selected features found in data columns for {virus}.")
                                continue
                            feature_cols = valid_feats
                            #print(f"    FS applied: {len(feature_cols)} features used.")

                        X_train = train_df[feature_cols]
                        y_train = train_df[target]
                        X_test = test_df[feature_cols]
                        y_test = test_df[target]
                        subject_ids = test_df["SubjectID"]

                        print(f"    Loaded Data ({virus}): Train Shape: {X_train.shape}, Test Shape: {X_test.shape}")
                        
                        if len(y_train.unique()) < 2:
                            continue
                            
                        # Train & Predict
                        model = get_model(classifier_name, seed)
                        try:
                            probs = train_predict(model, X_train, y_train, X_test)
                            
                            for subj, true_lbl, prob in zip(subject_ids, y_test, probs):
                                aggregated_predictions.append({
                                    "SubjectID": subj,
                                    "Virus": virus,
                                    "TrueLabel": true_lbl,
                                    "PredProb": prob
                                })
                        except Exception as e:
                            print(f"    Error training {virus}: {e}")
                    
                    # Save PP Files
                    if aggregated_predictions:
                        pp_df = pd.DataFrame(aggregated_predictions)
                        # Sort
                        try:
                            pp_df["_sort_key"] = pp_df["SubjectID"].apply(get_sort_key)
                            pp_df = pp_df.sort_values(by="_sort_key").drop(columns=["_sort_key"])
                        except Exception as e:
                            print(f"    Sorting failed: {e}")
                            
                        fs_suffix = f"_{fs_method}" if fs_method else ""
                        pp_filename = f"{t_time}_{classifier_name}_{level_name}_{target}{fs_suffix}_PP.csv"
                        pp_path = os.path.join(RESULTS_PP_DIR, pp_filename)
                        pp_df.to_csv(pp_path, index=False)
                        
                        # Metrics
                        try:
                            ap_score = average_precision_score(pp_df["TrueLabel"], pp_df["PredProb"])
                        except: ap_score = np.nan

                        try:
                            precision, recall, _ = precision_recall_curve(pp_df["TrueLabel"], pp_df["PredProb"])
                            auprc_score = auc(recall, precision)
                        except: auprc_score = np.nan
                        
                        preds = (pp_df["PredProb"] >= 0.5).astype(int)
                        acc = accuracy_score(pp_df["TrueLabel"], preds)
                        
                        # Confusion Matrix Metrics
                        try:
                            tn, fp, fn, tp = confusion_matrix(pp_df["TrueLabel"], preds, labels=[0, 1]).ravel()
                        except Exception as e:
                            print(f"    Error calculating CM: {e}")
                            tn, fp, fn, tp = 0, 0, 0, 0

                        try:
                            roc_auc = roc_auc_score(pp_df["TrueLabel"], pp_df["PredProb"])
                        except: roc_auc = np.nan
                        
                        print(f"    > Saved {pp_filename}. AUPRC: {auprc_score:.4f}, AP: {ap_score:.4f}, Acc: {acc:.4f}, TP: {tp}, FN: {fn}")
                        
                        summary_metrics.append({
                            "TimePoint": t_time,
                            "Classifier": classifier_name,
                            "Feature": level_name,
                            "Target": target,
                            "FS_Method": fs_method if fs_method else "None",
                            "AP": ap_score,
                            "AUPRC": auprc_score,
                            "Accuracy": acc,
                            "ROC_AUC": roc_auc,
                            "TP": tp,
                            "TN": tn,
                            "FP": fp,
                            "FN": fn,
                            "N_Samples": len(pp_df)
                        })

    # Save Summary
    if summary_metrics:
        summary_df = pd.DataFrame(summary_metrics)
        
        # Construct dynamic filename components
        c_str = classifiers[0] if len(classifiers) == 1 else "ALL"
        f_str = list(features.keys())[0] if len(features) == 1 else "ALL"
        t_str = targets[0] if len(targets) == 1 else "ALL"
        fs_str = fs_method if fs_method else "Baseline"
        
        out_name = f"Summary_Analysis_{c_str}_{f_str}_{t_str}_{fs_str}.csv"
        out_path = os.path.join(RESULTS_PP_DIR, out_name)
        summary_df.to_csv(out_path, index=False)
        print(f"\nSummary saved to {out_path}")
        
        try:
            print("\nGenerating Plots...")
            import plot_results
            plot_results.main() 
        except Exception as e:
            print(f"Plotting skipped/failed: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parameterized Viral Infection ML Analysis")
    
    parser.add_argument("--virus", type=str.upper, default="ALL",
                        help=f"Virus type: {', '.join(ALL_VIRUSES)} or ALL")
    
    parser.add_argument("--feature", type=str, default="ALL",
                        help=f"Feature type: {', '.join(ALL_LEVELS.keys())} or ALL (case sensitive)")
    
    parser.add_argument("--classifier", type=str, default="ALL",
                        help=f"Classifier: {', '.join(ALL_CLASSIFIERS)} or ALL")
    
    parser.add_argument("--target", type=str.upper, default="SC1",
                        choices=["SC1", "SC2", "ALL"],
                        help="Target variable: SC1 or SC2")
    
    parser.add_argument("--tp", type=str.upper, default="ALL",
                        help=f"Timepoint: {', '.join(map(str, TIME_POINTS))} or ALL")

    parser.add_argument("--fs", type=str, default=None,
                        help="Feature Selection Method (e.g., 'lasso'). Matches 'FS Method' column in CSV.")

    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    
    args = parser.parse_args()
    
    selected_viruses = ALL_VIRUSES if args.virus == "ALL" else [args.virus]
    
    if args.feature == "ALL":
        selected_features = ALL_LEVELS
    else:
        level_match = [k for k in ALL_LEVELS.keys() if k.lower() == args.feature.lower()]
        if level_match:
            selected_features = {level_match[0]: ALL_LEVELS[level_match[0]]}
        else:
            raise ValueError(f"Invalid feature: {args.feature}")

    selected_classifiers = ALL_CLASSIFIERS if args.classifier == "ALL" else [args.classifier]
    
    selected_targets = ["SC1", "SC2"] if args.target == "ALL" else [args.target]
    
    if args.tp == "ALL":
        selected_time_points = TIME_POINTS
    else:
        try:
            tp_val = int(args.tp)
            if tp_val in TIME_POINTS:
                selected_time_points = [tp_val]
            else:
                raise ValueError(f"Invalid Value: {tp_val}. Must be one of {TIME_POINTS}")
        except ValueError:
             raise ValueError(f"Invalid argument for --tp: {args.tp}")

    
    run_analysis(selected_classifiers, selected_targets, selected_viruses, selected_features, selected_time_points, args.seed, args.fs)
