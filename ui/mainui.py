from PyQt5 import QtWidgets, uic
from atoms import Atom, EnergyLevel, Transition
from atom_import import pickle_atom, load_from_pickle
from atom_library import populate_levels
from parsers import parse_NIST_levels
import sys
from warnings import warn

class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__()
        uic.loadUi('mainwindow.ui', self)

        self.loadedAtom: Atom
        # These should be dicts with "level", "F", and "m_F"
        self.ASSelectedLevel0: dict
        self.ASSelectedLevel1: dict
        # A list of the selected transitions? Will need to look into getting selected items from QTreeView and QListView
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
        if filename:
            if self.atomNameField.text() is (None or ""):
                # TODO: Make it so that the load button is just greyed out until there's a name specified
                warn("No atom name specified")
                return
            self.loadedAtom = Atom(name=self.atomNameField.text())
            populate_levels(df=parse_NIST_levels(filename), atom=self.loadedAtom,
                            I=self.doubleSpinBoxI.value(),
                            default_A=self.spinBoxDefaultA.value(),
                            default_B=self.spinBoxDefaultB.value(),
                            n_levels=self.spinBoxNumStates.value(),
                            hf_source=self.hfCSVField.text())
            self.loadedAtom.rezero()

    def hf_csv_file_browse(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self,
                                                            "Load a hyperfine csv", "", "CSV Files (*.csv)", options=options)
        # TODO: Make this actually load the hf data
        if filename:
            print(filename)

    def save_atom(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save atom", "", "Atom Files (*.atom)", options=options)
        if filename:
            pickle_atom(self.loadedAtom, filename=filename)

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
