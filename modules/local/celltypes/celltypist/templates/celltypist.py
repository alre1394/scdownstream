#!/usr/bin/env python3

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["CELLTYPIST_FOLDER"] = "./tmp/celltypist"

import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models as ct_models

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

models = "${models.join(' ')}".split()
save_probabilities = "${save_probabilities}" == "true"

adata_celltypist = adata.copy()  # make a copy of our adata
sc.pp.normalize_per_cell(
    adata_celltypist, counts_per_cell_after=10**4
)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

symbol_col = "${symbol_col}"
if symbol_col != "index" and symbol_col:
    if symbol_col not in adata_celltypist.var.columns:
        raise ValueError(f"Symbol column {symbol_col} not found in adata.var.columns")
    adata_celltypist.var_names = adata_celltypist.var[symbol_col]

# celltypist expects a string index, because it will make unique names by appending "-1", "-2"
# to duplicates if necessary. Cast other types (e.g. CategoricalIndex) to str:
adata_celltypist.var_names = adata_celltypist.var_names.astype(str)

df_list = []
probability_files = {}

# Define known organ atlas models (from https://www.celltypist.org/organs)
# Correct URL structure for all combined organ models: http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/ORGAN/models/Adult_Human_ORGAN.pkl
ORGAN_ATLAS_MODELS = {
    "Adult_Human_Blood": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Blood/models/Adult_Human_Blood.pkl",
    "Adult_Human_Bone_marrow": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Bone_marrow/models/Adult_Human_Bone_marrow.pkl",
    "Adult_Human_Heart": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Heart/models/Adult_Human_Heart.pkl",
    "Adult_Human_Hippocampus": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Hippocampus/models/Adult_Human_Hippocampus.pkl",
    "Adult_Human_Intestine": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Intestine/models/Adult_Human_Intestine.pkl",
    "Adult_Human_Kidney": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Kidney/models/Adult_Human_Kidney.pkl",
    "Adult_Human_Liver": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Liver/models/Adult_Human_Liver.pkl",
    "Adult_Human_Lung": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Lung/models/Adult_Human_Lung.pkl",
    "Adult_Human_Lymph_node": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Lymph_node/models/Adult_Human_Lymph_node.pkl",
    "Adult_Human_Pancreas": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Pancreas/models/Adult_Human_Pancreas.pkl",
    "Adult_Human_Skeletal_muscle": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Skeletal_muscle/models/Adult_Human_Skeletal_muscle.pkl",
    "Adult_Human_Spleen": "http://celltypist.cog.sanger.ac.uk/Resources/Organ_atlas/Spleen/models/Adult_Human_Spleen.pkl"
}

for model in models:
    print(f"Processing model: {model}")
    
    # Check if it's a file path (contains "/" or "\" or ends with .pkl)
    is_file_path = model.endswith(".pkl") or "/" in model or os.sep in model
    
    if is_file_path:
        # It's a local file path - load directly
        model_file = model
        model_name = os.path.basename(model).replace(".pkl", "")
        print(f"  Loading from file: {model_file}")
        print(f"  Current working directory: {os.getcwd()}")
        print(f"  Absolute path: {os.path.abspath(model_file)}")
        
        if not os.path.exists(model_file):
            # Check if it's an absolute path that needs to be made relative
            # Container work dir is isolated, so absolute paths won't work
            basename = os.path.basename(model_file)
            if os.path.exists(basename):
                model_file = basename
                print(f"  ✓ Found model in current directory: {basename}")
            else:
                print(f"  ✗ Model file not found in either location")
                print(f"    Tried: {model}")
                print(f"    Tried: {basename}")
                print(f"    Available files: {os.listdir('.')}")
                raise FileNotFoundError(
                    f"Model file not found: {model_file}. "
                    "When running in containers, absolute paths may not be accessible. "
                    "Use model names (for automatic download) or ensure the file is staged in the work directory."
                )
    
    elif model in ORGAN_ATLAS_MODELS:
        # It's an organ atlas model - download directly from organ atlas URL
        model_name = model
        model_file = f"{model_name}.pkl"
        organ_url = ORGAN_ATLAS_MODELS[model]
        
        print(f"  Downloading organ atlas model: {model_name}")
        print(f"  URL: {organ_url}")
        
        # Download the model file from organ atlas URL
        import urllib.request
        try:
            urllib.request.urlretrieve(organ_url, model_file)
            file_size_mb = os.path.getsize(model_file) / 1024**2
            print(f"  ✓ Organ atlas model downloaded successfully ({file_size_mb:.2f} MB)")
        except Exception as e:
            raise RuntimeError(
                f"Failed to download organ atlas model from {organ_url}: {e}. "
                "Please check your internet connection and the URL accessibility."
            )
    
    else:
        # It's a built-in celltypist model name - download using celltypist API
        model_name = model
        print(f"  Downloading built-in celltypist model: {model_name}")
        
        try:
            # Download the model (without .pkl extension - celltypist expects just the name)
            ct_models.download_models(model=model_name)
            model_file = f"{model_name}.pkl"
            print(f"  ✓ Built-in model downloaded successfully")
        except ValueError as e:
            # If model not found in built-in models, provide helpful error message
            print(f"  ✗ Model '{model_name}' not found in built-in celltypist models")
            print(f"  Available model types:")
            print(f"    - Built-in models: https://github.com/Teichlab/celltypist")
            print(f"    - Organ atlas models: {', '.join(ORGAN_ATLAS_MODELS.keys())}")
            print(f"    - Local file paths: /path/to/model.pkl")
            raise ValueError(f"Model '{model_name}' not recognized. {str(e)}")
    
    # Load the model
    print(f"  Loading model object from: {model_file}")
    model_obj = ct_models.Model.load(model_file)
    print(f"  ✓ Model loaded successfully")

    # Run celltypist annotation
    print(f"  Running celltypist annotation...")
    predictions = celltypist.annotate(
        adata_celltypist, model=model_obj
    )
    predictions_adata = predictions.to_adata()

    df_celltypist = predictions_adata.obs.loc[
        adata.obs.index, ["predicted_labels", "conf_score"]
    ]

    df_celltypist.columns = [f"celltypist:{model_name}", f"celltypist:{model_name}:conf"]
    df_list.append(df_celltypist)
    
    # Save full probability matrix if requested
    if save_probabilities:
        print(f"  Saving probability matrix to parquet file...")
        
        # Extract full probability matrix (already aligned with adata.obs.index)
        prob_matrix = predictions_adata.obsm["predicted_labels_probability"]
        
        # Get cell type names from the predictions
        cell_type_names = predictions_adata.obs['predicted_labels'].cat.categories.tolist()
        
        # Create DataFrame with cell barcodes as index and cell types as columns
        prob_df = pd.DataFrame(
            prob_matrix,
            index=adata.obs.index,
            columns=cell_type_names
        )
        
        # Save to compressed parquet file
        parquet_file = f"{prefix}_{model_name}_probabilities.parquet.gz"
        prob_df.to_parquet(parquet_file, compression='gzip', index=True)
        
        file_size_mb = os.path.getsize(parquet_file) / 1024**2
        print(f"  ✓ Probability matrix saved: {parquet_file} ({file_size_mb:.2f} MB)")
        
        # Store reference in adata metadata for later retrieval
        if "celltypist_probability_files" not in adata.uns:
            adata.uns["celltypist_probability_files"] = {}
        adata.uns["celltypist_probability_files"][model_name] = parquet_file
        probability_files[model_name] = parquet_file
    else:
        print(f"  Skipping probability matrix saving (celltypist_save_probabilities=False)")

df_celltypist = pd.concat(df_list, axis=1)
df_celltypist.to_pickle("${prefix}.pkl")

adata.obs = pd.concat([adata.obs, df_celltypist], axis=1)
adata.write_h5ad(f"{prefix}.h5ad")

# Save metadata about probability files
if save_probabilities and probability_files:
    metadata_df = pd.DataFrame({
        "model_name": list(probability_files.keys()),
        "parquet_file": list(probability_files.values())
    })
    metadata_df.to_csv(f"{prefix}_probabilities_metadata.csv", index=False)
    print(f"✓ Probability files metadata saved: {prefix}_probabilities_metadata.csv")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
        "celltypist": celltypist.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))