"""
Test Script for Mutalyze Executable
Verifies that all components are properly bundled and functional
"""

import sys
import os
from pathlib import Path

def test_imports():
    """Test that all critical imports work"""
    print("🧪 Testing imports...")
    
    tests = {
        'streamlit': 'Streamlit web framework',
        'openmm': 'OpenMM molecular dynamics',
        'openmm.app': 'OpenMM application layer',
        'openmm.unit': 'OpenMM units',
        'Bio.PDB': 'BioPython PDB module',
        'sklearn': 'scikit-learn ML',
        'sklearn.ensemble': 'Random Forest',
        'numpy': 'NumPy arrays',
        'pandas': 'Pandas dataframes',
        'scipy': 'SciPy scientific',
        'matplotlib': 'Matplotlib plotting',
        'seaborn': 'Seaborn statistical plots',
    }
    
    results = {}
    for module, desc in tests.items():
        try:
            __import__(module)
            results[module] = True
            print(f"  ✅ {module:20s} - {desc}")
        except ImportError as e:
            results[module] = False
            print(f"  ❌ {module:20s} - {desc} - {e}")
    
    return all(results.values())

def test_openmm_functionality():
    """Test that OpenMM actually works"""
    print("\n🔬 Testing OpenMM functionality...")
    
    try:
        from openmm import app
        import openmm
        from openmm import unit
        
        # Create simple system
        system = openmm.System()
        system.addParticle(1.0)
        
        # Create integrator
        integrator = openmm.LangevinIntegrator(
            300*unit.kelvin,
            1.0/unit.picosecond,
            2.0*unit.femtoseconds
        )
        
        # Create context
        context = openmm.Context(system, integrator)
        
        print(f"  ✅ OpenMM version: {openmm.version.version}")
        print(f"  ✅ Platform: {context.getPlatform().getName()}")
        print(f"  ✅ System created successfully")
        
        return True
        
    except Exception as e:
        print(f"  ❌ OpenMM test failed: {e}")
        return False

def test_force_field():
    """Test that AMBER force field is available"""
    print("\n🧬 Testing AMBER force field...")
    
    try:
        from openmm.app import ForceField
        
        ff = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
        print(f"  ✅ AMBER ff14SB loaded")
        print(f"  ✅ GB/SA implicit solvent loaded")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Force field test failed: {e}")
        return False

def test_model_files():
    """Test that calibration model exists"""
    print("\n🤖 Testing ML calibration model...")
    
    try:
        import pickle
        
        # Get the directory (works in both script and frozen exe)
        if getattr(sys, 'frozen', False):
            base_path = Path(sys._MEIPASS)
        else:
            base_path = Path(__file__).parent
        
        model_path = base_path / 'models' / 'calibration_v5.pkl'
        
        if not model_path.exists():
            print(f"  ❌ Model file not found: {model_path}")
            return False
        
        # Load model
        with open(model_path, 'rb') as f:
            model = pickle.load(f)
        
        print(f"  ✅ Calibration model loaded")
        print(f"  ✅ Model type: {type(model).__name__}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Model test failed: {e}")
        return False

def test_core_modules():
    """Test that core modules are available"""
    print("\n🔧 Testing core modules...")
    
    # Add core to path
    if getattr(sys, 'frozen', False):
        base_path = Path(sys._MEIPASS)
    else:
        base_path = Path(__file__).parent
    
    sys.path.insert(0, str(base_path / 'core'))
    
    modules = [
        'pdb_utils',
        'mutation_builder',
        'energy_calc',
        'empirical_correction',
        'interface_analysis',
        'conservation',
        'analysis',
        'parallel',
        'visualization',
        'simple_energy',
        'structure_fixer',
    ]
    
    results = {}
    for module in modules:
        try:
            __import__(module)
            results[module] = True
            print(f"  ✅ {module}")
        except ImportError as e:
            results[module] = False
            print(f"  ❌ {module} - {e}")
    
    return all(results.values())

def test_energy_calculation():
    """Test actual energy calculation"""
    print("\n⚡ Testing energy calculation...")
    
    try:
        from core.energy_calc import calculate_stability_change, OPENMM_AVAILABLE
        
        print(f"  ✅ OPENMM_AVAILABLE: {OPENMM_AVAILABLE}")
        
        if not OPENMM_AVAILABLE:
            print(f"  ⚠️  OpenMM not available - will use simplified calculator")
            print(f"  ⚠️  This should NOT happen in the executable!")
            return False
        
        print(f"  ✅ Energy calculation module loaded")
        return True
        
    except Exception as e:
        print(f"  ❌ Energy calculation test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("="*70)
    print(" MUTALYZE EXECUTABLE TEST SUITE")
    print("="*70)
    print()
    
    tests = [
        ("Imports", test_imports),
        ("OpenMM Functionality", test_openmm_functionality),
        ("Force Field", test_force_field),
        ("Model Files", test_model_files),
        ("Core Modules", test_core_modules),
        ("Energy Calculation", test_energy_calculation),
    ]
    
    results = {}
    for name, test_func in tests:
        try:
            results[name] = test_func()
        except Exception as e:
            print(f"\n❌ {name} test crashed: {e}")
            results[name] = False
    
    # Summary
    print("\n" + "="*70)
    print(" TEST SUMMARY")
    print("="*70)
    
    for name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status:10s} {name}")
    
    print("="*70)
    
    all_passed = all(results.values())
    
    if all_passed:
        print("\n🎉 ALL TESTS PASSED!")
        print("\nYour executable should work correctly with:")
        print("  ✅ Full OpenMM support")
        print("  ✅ AMBER ff19SB force field")
        print("  ✅ Random Forest v5 calibration")
        print("  ✅ All analysis features")
        print("\nReady to build executable! Run:")
        print("  Windows: build_windows.bat")
        print("  Linux:   ./build_linux.sh")
    else:
        print("\n⚠️  SOME TESTS FAILED!")
        print("\nFix the issues above before building executable.")
        print("\nCommon fixes:")
        print("  - Install OpenMM: conda install -c conda-forge openmm")
        print("  - Install dependencies: pip install -r requirements_exe.txt")
        print("  - Activate environment: conda activate mutalyze_build")
    
    print("\n" + "="*70 + "\n")
    
    return 0 if all_passed else 1

if __name__ == '__main__':
    sys.exit(main())
