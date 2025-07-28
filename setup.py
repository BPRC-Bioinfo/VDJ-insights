from setuptools import setup, find_packages

setup(
    name="vdj-insights",
    version="0.1.0",
    author="Jesse Mittertreiner, Sayed Jamiel Mohammadi, Giang Le, Jesse Bruijnesteijn, Suzan Ott",
    author_email="jaimymohammadi@gmail.com",
    description="VDJ-Insights provides a robust framework for the accurate annotation of complex genomic immune regions.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(include=["vdj_insights", "vdj_insights.*"]),
    include_package_data=True,
    package_data={"vdj_insights": ["**/*"]},
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires=">=3.7",
    install_requires=[
        "biopython==1.85",
        "openpyxl==3.1.5",
        "PyYAML",
        "tqdm==4.67.1",
        "psutil==7.0.0",
        "matplotlib==3.10.3",
        "seaborn==0.13.2",
        "bs4==0.0.2",
        "venny4py==1.0.3",
        "requests==2.32.4",
        "Flask==3.1.1",
        "flask_caching==2.3.1",
        "bokeh==3.7.3",
        "plotly==6.2.0",
        "matplotlib-venn==1.1.2",
        "dna-features-viewer==3.1.5",
        "ghostscript==0.8.1"
    ],
    entry_points={
        "console_scripts": [
            "vdj-insights=vdj_insights.scripts.vdj_insights:main",
        ],
    },
)
