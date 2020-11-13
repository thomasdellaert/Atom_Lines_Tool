import os.path
import pickle
# import jsons

# Currently, json saving doesn't work

# def dump_to_json(atom, filename=None):
#     if filename is None:
#         filename = atom.name
#     json_string = jsons.dump(atom, verbose=jsons.Verbosity.WITH_EVERYTHING)
#     file = open(filename+".json", "w")
#     file.write(str(json_string))
#     file.close()
#
# def load_from_json(file):
#     with open(file) as f:
#         data = f.read()
#     return jsons.load(data)

def pickle_atom(atom, filename=None):
    if filename is None:
        filename = atom.name
    file = open(filename+".atom", "wb")
    pickle.dump(atom, file)
    file.close()

def load_from_pickle(filename):
    file = open(filename, "rb")
    p = pickle.load(file)
    file.close()
    return p


if __name__ == "__main__":
    print(True)
