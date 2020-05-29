from setuptools import setup


setup(
    name = "NICD",
    packages = ["NICD"],
    install_requires = ["numpy", "scipy", "tqdm", "networkx", "pandas", "ndex2"],
    version = "0.0.1",
    description = "Identifying potential common gene drivers that can cause multiple-diseases.",
    author="Hyun-Hwan Jeong",
    author_email="hyun-hwan.jeong@bcm.edu",
    keywords=["disease", "gene", "regulation", "bioinformatics"],
    url="https://github.com/hyunhwaj/NICD",
    scripts=[
        "bin/NICD",
        "bin/NICD_DB_download"
    ],
    include_package_data=True
)
