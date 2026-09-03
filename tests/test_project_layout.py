"""Repository-distribution contract tests."""

from pathlib import Path


REPO = Path(__file__).resolve().parent.parent


def test_repository_is_distributed_as_a_snakemake_workflow(pytestconfig):
    """Do not advertise a broken Python wheel for a repository-run workflow."""

    assert (REPO / "workflow" / "Snakefile").is_file()
    assert pytestconfig.inipath == REPO / "pytest.ini"
    assert (REPO / "ruff.toml").is_file()
    assert not any(
        (REPO / filename).exists()
        for filename in ("setup.py", "setup.cfg", "pyproject.toml")
    )
