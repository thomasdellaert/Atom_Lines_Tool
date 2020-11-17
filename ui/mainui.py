from PyQt5 import QtWidgets, uic
from atoms import Atom, EnergyLevel, Transition
from atom_import import pickle_atom, load_from_pickle
import sys


class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__()
        uic.loadUi('mainwindow.ui', self)

        self.loadedAtom: Atom
        # These should be dicts with "level", "F", and "m_F"
        self.ASSelectedLevel0: dict
        self.ASSelectedLevel1: dict
        self.ASSelectedTransitions: list

        # Define the Atom Setup tab buttons

        self.loadNISTCSVButton.clicked.connect(self.nist_csv_file_browse)
        self.browseCSVButton.clicked.connect(self.hf_csv_file_browse)
        self.saveButton.clicked.connect(self.save_atom)
        self.createTransitionButton.clicked.connect(self.create_transition)
        # self.deleteTransitionButton.clicked.connect(self.delete_transition)

        # Define the

        self.show()

    def nist_csv_file_browse(self):
        """opens a file browser from which to load a NIST csv"""
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load a NIST csv", "", "CSV Files (*.csv)", options=options)
        # TODO: Make this actually load the atom
        if filename:
            # self.
            print(filename)

    def hf_csv_file_browse(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self,
                                                            "Load a hyperfine csv", "", "CSV Files (*.csv)", options=options)
        # TODO: Make this actually load the hf data
        if filename:
            print(filename)

    def save_atom(self, loaded_atom):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save atom", "", "Atom Files (*.atom)", options=options)
        # TODO: Make this actually pickle the atom as loaded

        if filename:
            print(filename)

    def create_transition(self):
        try:
            t = Transition(level_0=self.ASSelectedLevel0["level"], F_0=self.ASSelectefLevel0["F"], m_F_0=self.ASSelectedLevel0["m_F"],
                           level_1=self.ASSelectedLevel1["level"], F_1=self.ASSelectefLevel1["F"], m_F_1=self.ASSelectedLevel1["m_F"])
            self.loadedAtom.add_transition(t)
        except:
            pass
        # TODO: Make this work
        print("creating transition")

    def delete_transition(self):
        self.loadedAtom.remove_transition(self.ASSelectedTransitions, mode="transitions")
        # TODO: Add this button and make it work


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = Ui()
    app.exec()
