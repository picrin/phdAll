import nbformat

from nbconvert.preprocessors import ExecutePreprocessor

def executeNotebook(notebookName, wd):
    with open(notebookName) as file:
        nb = nbformat.read(file, as_version=4)
        ep = ExecutePreprocessor(timeout=600, kernel_name='ir')
        ep.preprocess(nb, {'metadata': {'path': wd}})
