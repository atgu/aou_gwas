from setuptools import find_packages, setup

setup(
    name="aou_gwas",
    version="0.0.1",
    packages=find_packages(),
        install_requires=[
        "google-cloud-storage",
        "gnomad",
        "hail",
        "numpy",
        "pandas",
        "tqdm",
    ],
)
