from setuptools import setup, find_packages

setup(
    name="bctool",
    version="1.0.0",
    description="Bisulfite conversion alignment & benchmarking toolkit",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "click>=8.0",
        "pyyaml>=5.1",
        "pandas>=1.3",
        "numpy>=1.21",
        "biopython>=1.79",
        "matplotlib>=3.4",
        "seaborn>=0.11",
        "rich>=10.0",
    ],
    entry_points={
        "console_scripts": [
            "bctool=bctool.cli:main",
        ],
    },
    python_requires=">=3.8",
)
