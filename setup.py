from setuptools import setup, find_packages

setup(
    name="vdj-insights",
    version="0.1.111",
    author="Jesse Mittertreiner, Sayed Jamiel Mohammadi, Giang Le",
    author_email="jaimymohammadi@gmail.com",
    description="VDJ insights offers a robust framework for analyzing, assembling, and annotating long sequence reads from Pacific Biosciences (PacBio) and Oxford Nanopore Technologies (ONT).",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(include=["vdj_insights", "vdj_insights.*"]),
    include_package_data=True,
    package_data={"vdj_insights": ["**/*"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        "biopython",
        "PyYAML",
        "tqdm",
        "psutil",
        "matplotlib",
        "seaborn",
        "bs4",
        "venny4py"
    ],
    entry_points={
        "console_scripts": [
            "vdj-insights=vdj_insights.scripts.vdj_insights:main",
        ],
    },
)
