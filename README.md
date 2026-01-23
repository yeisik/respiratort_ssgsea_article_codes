# GSE73072 Dataset Construction (RMA + ComBat + Virus-wise Train/Test Splits)

This page documents how we reconstructed the datasets used in our paper from the public GEO study **GSE73072** (Affymetrix Human Genome U133A 2.0). The goal is to provide a **fully reproducible** pipeline: starting from raw `.CEL` files and ending with **virus-specific, timepoint-specific train/test CSV files** used for infection and symptom prediction.

---

## 📋 Overview of the Pipeline

We follow three main steps:

1.  **RMA normalization from raw CEL files**
    * Implemented in: `celread.r`
2.  **Batch effect correction (ComBat)** across the 7 experiments
    * Implemented in: `combat.ipynb`
3.  **Train/Test splitting per experiment and virus-wise dataset assembly**
    * Outputs saved under: `Datasets/`
    * Labels stored in: `final_all_labels.csv`
4.  **Machine Learning Analysis Script (
run_ml_analysis_v2.py
)
**

---

## Step 1: RMA Normalization from Raw CEL Files

* **Input:** `GSE73072 RAW/` folder containing Affymetrix `.CEL` files
* **Script:** `celread.r`
* **Outputs:** Probe-level and/or gene-level expression matrices (depending on the selected CDF mapping)

Affymetrix `.CEL` files contain **raw probe intensities**, which are not directly comparable across arrays due to technical effects (scanner behavior, background noise, labeling differences, etc.). Therefore, we apply **Robust Multi-array Average (RMA)**, a standard preprocessing pipeline that improves comparability and reduces noise by performing:

* **Background correction**
* **Quantile normalization** (forces comparable intensity distributions across arrays)
* **Summarization** (combines multiple probes into a single expression value per probe set / gene)

RMA yields a stable expression matrix that is suitable for downstream machine learning and is widely used in the microarray literature.

> **Note:** In our workflow, RMA is the first essential step to remove within-array technical variation and to obtain consistent expression values from raw CEL intensities.

---

## Step 2: Batch Effect Correction Across Experiments

* **Input:** RMA-normalized expression matrix/matrices from Step 1
* **Notebook:** `combat.ipynb`
* **Outputs:** ComBat-corrected expression matrices

GSE73072 is composed of **seven distinct experiments** (DEE1–DEE5, HRV UVA, HRV DUKE). Even after RMA, aggregated datasets can still exhibit **between-experiment batch effects** (systematic shifts caused by different experimental dates, personnel, protocols, and handling).

To reduce these cross-experiment artifacts, we apply **ComBat** (empirical Bayes batch correction) using the experiment identifier as the **batch variable**. This step harmonizes expression distributions across experiments while preserving biological signal relevant to infection/symptom outcomes.

> **Summary:** **RMA** addresses within-array/within-batch normalization, while **ComBat** targets *between-experiment* technical variation introduced when multiple experiments are pooled.

---

## Step 3: Train/Test Splits and Virus-wise Dataset Assembly

After obtaining ComBat-corrected expression matrices, we construct the supervised learning datasets used in the paper.

### 3.1 Experiments and Virus Mapping

GSE73072 consists of the following **7 experiments** corresponding to **4 respiratory viruses**:

| Experiment | Virus Type |
| :--- | :--- |
| **RSV DEE1** | RSV |
| **H3N2 DEE2** | H3N2 |
| **H1N1 DEE3** | H1N1 |
| **H1N1 DEE4** | H1N1 |
| **H3N2 DEE5** | H3N2 |
| **HRV UVA** | HRV |
| **HRV DUKE** | HRV |

### 3.2 Train/Test Splitting

We first split each experiment into **train** and **test** sets (as described in the paper), then merge experiments that belong to the same virus to form **virus-wise datasets**.

This yields one train/test dataset per virus (instead of an experiment-wise analysis), enabling consistent virus-specific model evaluation.

### 3.3 Output Datasets Directory Structure

All finalized datasets are stored under `Datasets/` using the following structure:

```text
Datasets/
 └── VirusName/
      └── FeatureRepresentationName/
           ├── VirusName_train_Timepoint.csv
           └── VirusName_test_Timepoint.csv
#### Example File Names
As found in this repository:
* `H1N1_Probe_test_T0.csv`
* `H1N1_Probe_train_T24.csv`

#### Nomenclature Key
* **VirusName:** `{H1N1, H3N2, HRV, RSV}`
* **FeatureRepresentationName:** Corresponds to the representation used in the paper (e.g., `Probe`, `Gene`, or others).
* **Timepoint:** Standardized evaluation category (e.g., `T0`, `T24`, `T48`, `T72`, `T96`, `T120`).

### 📂 Labels File
All class labels used in the paper are provided in the root directory:

`final_all_labels.csv`

This file contains the ground-truth labels required for:
* **Infection prediction** (infected vs non-infected)
* **Symptom prediction** (symptomatic vs non-symptomatic; according to the paper’s definition)
```

### 4. Machine Learning Analysis Script (`run_ml_analysis_v2.py`)

This script performs machine learning analysis on viral infection data. It supports various classifiers, feature sets, and optional feature selection filtering.

## Usage

Run the script from the command line using the following parameters:

```bash
python run_ml_analysis_v2.py [OPTIONS]


### Parameters

| Argument | Description | Options | Default |
| :--- | :--- | :--- | :--- |
| `--virus` | Specify the virus type to analyze. | `H1N1`, `H3N2`, `HRV`, `RSV`, `ALL` | `ALL` |
| `--feature` | Specify the feature set (representation) to use. | `Probe`, `ssGSEA`, `Combined`, `ALL` | `ALL` |
| `--classifier` | Specify the machine learning classifier. | `LR`, `XGB`, `SVM`, `KNN`, `NuSVC`, `GNB`, `RF`, `ALL` | `ALL` |
| `--target` | Specify the target variable (Label). | `SC1` (Infection), `SC2` (Symptoms), `ALL` | `SC1` |
| `--tp` | Specify the timepoint (hours post-inoculation). | `0`, `24`, `48`, `72`, `96`, `120`, `ALL` | `ALL` |
| `--fs` | Apply Feature Selection filtering. Must match a method name in your FS results CSV. | Any string (e.g., `lasso`, `RF*`, `ttest`) | `None` (Baseline) |
| `--seed` | Set the random seed for reproducibility. | Integer | `42` |

### Examples

**1. Run Baseline Analysis (No Feature Selection)**
Run Random Forest on H1N1 data using Combined features for all timepoints, predicting Infection (SC1).
```bash
python run_ml_analysis_v2.py --virus H1N1 --feature Combined --classifier RF --target SC1
```

**2. Run Feature Selection Analysis**
Apply 'lasso' feature selection. The script will look up the selected features for each combination in `filtered_FS_results.csv` and filter the data accordingly.
```bash
python run_ml_analysis_v2.py --virus H1N1 --feature Combined --classifier RF --target SC1 --fs lasso
```

**3. Run Specific Timepoint**
Analyze only T96 with Support Vector Machine.
```bash
python run_ml_analysis_v2.py --tp 96 --classifier SVM --target SC2
```

**4. Run All Configurations**
Run all viruses, all features, all classifiers, and all timepoints.
```bash
python run_ml_analysis_v2.py --virus ALL --feature ALL --classifier ALL --target ALL --tp ALL
```

### Output
Results are saved in the `Results_PP` directory:
- **Prediction Probabilities**: `[Time]_[Classifier]_[Feature]_[Target]_[FSMethod]_PP.csv`
- **Summary Metrics**: `Summary_Analysis_[Classifier]_[Feature]_[Target]_[FSMethod].csv` containing Accuracy, AUPRC, AP, ROC AUC, and Confusion Matrix metrics (TP, TN, FP, FN).
