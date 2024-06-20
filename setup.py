import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SimulatedTumorData",
    version="0.0.1",
    author="Kristy Schlueter-Kuck and Claudia Chu",
    author_email="cchu@broadinstitute.org",
    description="Package to generate simulated tumor data and initial example data",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/getzlab/SimulatedTumorData",
    project_urls={
        "Bug Tracker": "https://github.com/getzlab/SimulatedTumorData/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "."},
    packages=['src'], #setuptools.find_packages(where="."),
    python_requires=">=3.6",
    install_requires = [
        'AnnoMate>=0.0.2',
        'seaborn',
        'dash-cytoscape',
        'dash',
        'pysam',
        'scipy',
        'wheel',
    ]
)   

