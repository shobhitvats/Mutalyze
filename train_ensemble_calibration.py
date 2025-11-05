#!/usr/bin/env python3
"""
Train improved calibration using ensemble methods and cross-validation
Optimizes for oracle-level performance (r>0.8, RMSE<0.5)
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler
import pickle
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def load_training_data(log_file='training_v5_with_fixer.log'):
    """Extract successful predictions from training log"""
    
    predictions = []
    
    with open(log_file, 'r') as f:
        for line in f:
            if '✓ Exp:' in line and 'Raw:' in line:
                try:
                    parts = line.split('Exp:')[1].split('Raw:')
                    exp = float(parts[0].strip().rstrip(','))
                    raw = float(parts[1].strip())
                    
                    # Only include finite values
                    if np.isfinite(raw) and np.isfinite(exp):
                        predictions.append({'experimental': exp, 'raw': raw})
                except:
                    continue
    
    df = pd.DataFrame(predictions)
    logger.info(f"Loaded {len(df)} successful predictions with finite values")
    return df


def create_features(raw_values):
    """Create features for non-linear calibration"""
    X = np.array(raw_values).reshape(-1, 1)
    
    # Add polynomial features
    X_poly = np.column_stack([
        X,                    # Linear
        X**2,                 # Quadratic
        X**3,                 # Cubic
        np.abs(X),            # Absolute value
        np.sign(X) * np.sqrt(np.abs(X)),  # Square root with sign
    ])
    
    return X_poly


def train_ensemble_model(df):
    """Train ensemble calibration model with cross-validation"""
    
    X = create_features(df['raw'].values)
    y = df['experimental'].values
    
    logger.info(f"\nTraining on {len(df)} mutations")
    logger.info(f"Raw ΔΔG range: [{df['raw'].min():.1f}, {df['raw'].max():.1f}]")
    logger.info(f"Experimental ΔΔG range: [{df['experimental'].min():.1f}, {df['experimental'].max():.1f}]")
    
    # Test multiple models
    models = {
        'Ridge': Ridge(alpha=10.0),
        'RandomForest': RandomForestRegressor(
            n_estimators=100,
            max_depth=5,
            min_samples_split=5,
            min_samples_leaf=2,
            random_state=42
        ),
        'GradientBoosting': GradientBoostingRegressor(
            n_estimators=50,
            max_depth=3,
            learning_rate=0.1,
            random_state=42
        )
    }
    
    results = {}
    
    for name, model in models.items():
        logger.info(f"\n{'='*60}")
        logger.info(f"Training {name}")
        logger.info(f"{'='*60}")
        
        # Cross-validation
        cv = KFold(n_splits=5, shuffle=True, random_state=42)
        cv_scores = cross_val_score(model, X, y, cv=cv, scoring='r2')
        cv_r2 = cv_scores.mean()
        
        # Train on full data
        model.fit(X, y)
        y_pred = model.predict(X)
        
        # Metrics
        r, _ = pearsonr(y, y_pred)
        rmse = np.sqrt(np.mean((y - y_pred)**2))
        mae = np.mean(np.abs(y - y_pred))
        
        logger.info(f"Training R²: {model.score(X, y):.4f}")
        logger.info(f"Training r: {r:.4f}")
        logger.info(f"Training RMSE: {rmse:.2f} kcal/mol")
        logger.info(f"Training MAE: {mae:.2f} kcal/mol")
        logger.info(f"CV R² (5-fold): {cv_r2:.4f} ± {cv_scores.std():.4f}")
        
        results[name] = {
            'model': model,
            'cv_r2': cv_r2,
            'train_r': r,
            'train_rmse': rmse,
            'train_mae': mae
        }
    
    # Select best model based on CV R²
    best_name = max(results.keys(), key=lambda k: results[k]['cv_r2'])
    best_model = results[best_name]['model']
    
    logger.info(f"\n{'='*60}")
    logger.info(f"BEST MODEL: {best_name}")
    logger.info(f"{'='*60}")
    logger.info(f"CV R²: {results[best_name]['cv_r2']:.4f}")
    logger.info(f"Training r: {results[best_name]['train_r']:.4f}")
    logger.info(f"Training RMSE: {results[best_name]['train_rmse']:.2f} kcal/mol")
    
    return best_model, best_name, results


def save_model(model, model_name, filename='models/calibration_v5.pkl'):
    """Save trained model"""
    import os
    os.makedirs('models', exist_ok=True)
    
    model_data = {
        'model': model,
        'model_name': model_name,
        'version': 'v5.0',
        'description': 'Ensemble calibration with polynomial features'
    }
    
    with open(filename, 'wb') as f:
        pickle.dump(model_data, f)
    
    logger.info(f"\nModel saved to {filename}")


def generate_calibration_function(model, model_name):
    """Generate Python function for the calibration"""
    
    logger.info(f"\n{'='*60}")
    logger.info("CALIBRATION FUNCTION")
    logger.info(f"{'='*60}")
    
    if model_name == 'Ridge':
        # Extract coefficients
        coef = model.coef_
        intercept = model.intercept_
        
        logger.info(f"""
def apply_correction_v5(ddg_raw):
    \"\"\"
    Apply v5 polynomial calibration (Ridge regression)
    
    Formula:
        ΔΔG = {intercept:.6f} + 
              {coef[0]:.6f} * raw +
              {coef[1]:.6f} * raw² +
              {coef[2]:.6f} * raw³ +
              {coef[3]:.6f} * |raw| +
              {coef[4]:.6f} * sign(raw) * √|raw|
    \"\"\"
    import numpy as np
    
    raw = ddg_raw
    features = [
        raw,
        raw**2,
        raw**3,
        np.abs(raw),
        np.sign(raw) * np.sqrt(np.abs(raw))
    ]
    
    return {intercept:.6f} + np.dot([{', '.join(f'{c:.6f}' for c in coef)}], features)
""")
    
    else:
        logger.info(f"""
def apply_correction_v5(ddg_raw):
    \"\"\"
    Apply v5 calibration using {model_name}
    
    Note: Requires loading pickled model from models/calibration_v5.pkl
    \"\"\"
    import pickle
    import numpy as np
    
    # Load model
    with open('models/calibration_v5.pkl', 'rb') as f:
        model_data = pickle.load(f)
    
    model = model_data['model']
    
    # Create features
    raw = ddg_raw
    X = np.array([[
        raw,
        raw**2,
        raw**3,
        np.abs(raw),
        np.sign(raw) * np.sqrt(np.abs(raw))
    ]])
    
    return model.predict(X)[0]
""")


if __name__ == '__main__':
    logger.info("="*60)
    logger.info("ENSEMBLE CALIBRATION TRAINING (v5)")
    logger.info("="*60)
    
    # Load data
    df = load_training_data()
    
    if len(df) < 10:
        logger.error("Insufficient training data (need at least 10 samples)")
        exit(1)
    
    # Train ensemble model
    model, model_name, results = train_ensemble_model(df)
    
    # Save model
    save_model(model, model_name)
    
    # Generate calibration function
    generate_calibration_function(model, model_name)
    
    logger.info("\n" + "="*60)
    logger.info("TRAINING COMPLETE")
    logger.info("="*60)
    logger.info(f"Best model: {model_name}")
    logger.info(f"Ready to deploy v5 calibration")
