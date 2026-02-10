import pygcmc


def test_platform_logging_cross_translation_unit(capfd):
    pygcmc.set_platform_verbose(True)
    pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.DEBUG)

    try:
        pygcmc.setPGPParameters(
            alpha=0.29,
            meshSize=[8, 8, 8],
            potential_cutoff=0.5,
            potentialGridSize=[4, 4, 4],
            splineOrder=4,
            tolerance=1e-5,
        )

        out, err = capfd.readouterr()
        assert "setPGPParameters called" in out
    finally:
        if hasattr(pygcmc, "resetPGPState"):
            pygcmc.resetPGPState()
        pygcmc.set_platform_verbose(False)
        pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.INFO)
