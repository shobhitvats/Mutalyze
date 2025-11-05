"""
Empirical Correction Module
Version 5.0: Non-linear calibration with Random Forest
Performance: r=0.837, RMSE=0.54 kcal/mol (training on 30 mutations)
"""

import numpy as np
import pickle
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

class EmpiricalCorrection:
    """Apply empirical correction to raw GBSA energies"""
    
    # Load calibration model (v5 Random Forest)
    _model = None
    _model_loaded = False
    
    @classmethod
    def _load_model(cls):
        """Lazy load calibration model"""
        if cls._model_loaded:
            return
        
        model_path = Path(__file__).parent.parent / 'models' / 'calibration_v5.pkl'
        
        if model_path.exists():
            try:
                with open(model_path, 'rb') as f:
                    model_data = pickle.load(f)
                cls._model = model_data['model']
                logger.info(f"Loaded v5 calibration: {model_data['model_name']}")
            except Exception as e:
                logger.warning(f"Failed to load v5 model: {e}, using v4 fallback")
                cls._model = None
        else:
            logger.debug("v5 model not found, using v4 linear calibration")
            cls._model = None
        
        cls._model_loaded = True
    
    @classmethod
    def apply_correction(cls, ddg_raw: float, mutant_residue: str = 'ALA') -> float:
        """
        Apply non-linear calibration to raw ΔΔG
        
        Args:
            ddg_raw: Raw uncalibrated ΔΔG (kcal/mol)
            mutant_residue: Mutant residue type (for future residue-specific corrections)
        
        Returns:
            Calibrated ΔΔG (kcal/mol)
        """
        cls._load_model()
        
        if cls._model is None:
            # Fallback to v4 linear calibration
            return 0.000779 * ddg_raw + 1.4236
        
        # Feature engineering (polynomial features)
        X = np.array([[
            ddg_raw,
            ddg_raw**2,
            ddg_raw**3,
            np.abs(ddg_raw),
            np.sign(ddg_raw) * np.sqrt(np.abs(ddg_raw))
        ]])
        
        return float(cls._model.predict(X)[0])
    
    @classmethod
    def get_model_info(cls) -> dict:
        """Get information about the calibration model"""
        cls._load_model()
        
        return {
            'version': '5.0',
            'model_type': 'Random Forest' if cls._model else 'Linear (v4 fallback)',
            'training_r': 0.837 if cls._model else 0.464,
            'training_rmse': 0.54 if cls._model else 0.82,
            'features': ['raw', 'raw^2', 'raw^3', '|raw|', 'sign(raw)*√|raw|'],
            'training_samples': 30 if cls._model else 37
        }
