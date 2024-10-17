import dill

with open('classifier.pkl', 'rb') as pkl:
    test = dill.load(pkl)

    dill.source.getsource(test)