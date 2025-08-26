# tests/simulation/movementCPP/movement_module_basic_funcs.py
"""Movement module basic functionality tests - extracted functions."""

import pytest
import pygcmc

@pytest.fixture
def setup_system():
    """Setup a basic test system."""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.5]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    
    return state, params

def test_module_creation_basic():
        """Test MovementModule creation with different constructors."""
        # Default constructor
        mover1 = pygcmc.movement.MovementModule()
        assert mover1 is not None
        
        # Constructor with params
        params = pygcmc.movement.MovementParams()
        mover2 = pygcmc.movement.MovementModule()
        mover2.setParams(params)
        assert mover2 is not None
    
def test_parameter_setting_and_getting():
        """Test setting and getting parameters."""
        mover = pygcmc.movement.MovementModule()
        
        # Set parameters
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.chemicalPotential = -20.0
        mover.setParams(params)
        
        # Get parameters
        params_back = mover.getParams()
        assert params_back.temperature == pytest.approx(300.0)
        assert params_back.chemicalPotential == pytest.approx(-20.0)
    
def test_insertion_basic(setup_system):
        """Test basic insertion functionality."""
        state, params = setup_system
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Attempt insertion
        result = mover.attemptInsertion(state)
        
        # Check result
        assert isinstance(result.accepted, bool)
        assert result.moveType == "insert"  # Insertion
        assert isinstance(result.energyChange, float)
    
def test_deletion_basic(setup_system):
        """Test basic deletion functionality."""
        state, params = setup_system
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # First insert some atoms
        inserted = False
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Try deletion
            result = mover.attemptDeletion(state)
            assert isinstance(result.accepted, bool)
            assert result.moveType == "delete"  # Deletion
    
def test_translation_basic(setup_system):
        """Test basic translation functionality."""
        state, params = setup_system
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # First insert an atom
        inserted = False
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Try translation
            result = mover.attemptTranslation(state)
            assert isinstance(result.accepted, bool)
            assert result.moveType == "translate"  # Translation
    
def test_rotation_basic(setup_system):
        """Test basic rotation functionality."""
        state, params = setup_system
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # First insert an atom
        inserted = False
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Try rotation
            result = mover.attemptRotation(state)
            assert isinstance(result.accepted, bool)
            assert result.moveType == "rotate"  # Rotation
    
def test_cavity_bias_insertion(setup_system):
        """Test cavity-biased insertion."""
        state, params = setup_system
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        params.useCavityBias = True
        mover.setParams(params)
        
        # Attempt cavity-biased insertion
        result = mover.attemptCavityBiasInsertion(state)
        
        assert isinstance(result.accepted, bool)
        assert result.moveType == "insert"  # Still insertion
    
def test_config_bias_rotation_basic(setup_system):
        """Test configurational bias rotation basic functionality."""
        state, params = setup_system
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        params.useConfigBias = True
        mover.setParams(params)
        
        # First insert an atom
        inserted = False
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Try config bias rotation
            result = mover.attemptConfigBiasRotation(state)
            assert isinstance(result.accepted, bool)
    
def test_acceptance_rate_calculation(setup_system):
        """Test acceptance rate calculation."""
        state, params = setup_system
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Reset statistics
        mover.resetStatistics()
        
        # Perform moves and track
        accepts = 0
        attempts = 100
        
        for _ in range(attempts):
            result = mover.attemptInsertion(state)
            if result.accepted:
                accepts += 1
        
        # Calculate acceptance rate
        rate = mover.calculateAcceptanceRate("insert")
        
        # Should be close to manual calculation
        expected_rate = accepts / attempts
        assert abs(rate - expected_rate) < 0.1
    
def test_find_cavities_integration(setup_system):
        """Test cavity finding integration."""
        state, params = setup_system
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        params.useCavityBias = True
        mover.setParams(params)
        
        # Find cavities
        cavities = mover.findCavities(state)
        
        # Should return a list of positions
        assert isinstance(cavities, list)
        if len(cavities) > 0:
            # Check first cavity
            cavity = cavities[0]
            assert hasattr(cavity, 'x')
            assert hasattr(cavity, 'y')
            assert hasattr(cavity, 'z')
            
            # Should be within box
