from distutils.core import setup

setup(
    name = "NICD",
    packages = ["ncid"],
    install_requires = ["numpy", "scipy", "tqdm", "sklearn", "pandas"],
    version = "0.0.1",
    description = "Identifying potential common gene drivers that can cause multiple-diseases.",
    author="Hyun-Hwan Jeong",
    author_email="hyun-hwan.jeong@bcm.edu",
    keywords=["disease", "gene", "regulation", "bioinformatics"],
    url="https://github.com/hyunhwaj/NICD",
    scripts=[
        "bin/NICD"
    ]
)