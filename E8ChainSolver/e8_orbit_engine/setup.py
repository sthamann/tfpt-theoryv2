#!/usr/bin/env python
"""Setup script for E8 Orbit Engine."""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text() if readme_path.exists() else ""

setup(
    name="e8-orbit-engine",
    version="0.1.0",
    author="Stefan Hamann",
    author_email="stefan.hamann@example.com",
    description="E8 Orbit Engine - Quadratic damping from nilpotent orbits",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/stefanhamann/e8-orbit-engine",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.24.0",
        "pandas>=2.0.0",
        "scipy>=1.11.0",
        "mpmath>=1.3.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "typer>=0.9.0",
        "jinja2>=3.1.0",
        "openpyxl>=3.1.0",
        "rich>=13.0.0",
        "plotly>=5.17.0",
        "kaleido>=0.2.1",
    ],
    extras_require={
        "dev": [
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",
            "black>=23.7.0",
            "ruff>=0.0.287",
            "mypy>=1.5.0",
            "ipykernel>=6.25.0",
            "jupyterlab>=4.0.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "e8-orbit=e8_orbit_engine.cli:app",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords="E8 orbit nilpotent physics damping",
)
