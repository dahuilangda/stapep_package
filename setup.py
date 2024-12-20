from setuptools import setup

setup(
    name="stapep",
    version="0.2",
    description="A toolset for stapled peptides.",
    author="dahuilangda",
    author_email="dahuilangda@hotmail.com",
    packages=["stapep", 
              "stapep.models", 
              "stapep.templates.AIB",
              "stapep.templates.NLE",
              "stapep.templates.PR3",
              "stapep.templates.PR5",
              "stapep.templates.PR8",
              "stapep.templates.PS3",
              "stapep.templates.PS5",
              "stapep.templates.PS8",],
    package_data={
        "stapep.models": ["lgb_model.sav"],
        "stapep.templates.AIB": ['*'],
        "stapep.templates.NLE": ['*'],
        "stapep.templates.PR3": ['*'],
        "stapep.templates.PR5": ['*'],
        "stapep.templates.PR8": ['*'],
        "stapep.templates.PS3": ['*'],
        "stapep.templates.PS5": ['*'],
        "stapep.templates.PS8": ['*'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)