from PyQt5 import QtWidgets, uic, QtCore
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
        self.loadAtomButton.clicked.connect(self.atom_file_browse)
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

    def atom_file_browse(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self,
                                                            "Load an atom", "", "Atom Files (*.atom)", options=options)
        if filename:
            self.loadedAtom = load_from_pickle(filename=filename)

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


# TODO Implement a model of the atom for views to act on
class AtomModel(QtCore.QAbstractItemModel):
    def __init__(self, atom, *args, **kwargs):
        super(AtomModel, self).__init__(*args, **kwargs)
        self.atom = atom
        self._root = CustomNode(None)

    def addChild(self, in_node, in_parent):
        if not in_parent or not in_parent.isValid():
            parent = self._root
        else:
            parent = in_parent.internalPointer()
        parent.addChild(in_node)

    def index(self, in_row, in_column, in_parent=None):
        if not in_parent or not in_parent.isValid():
            parent = self._root
        else:
            parent = in_parent.internalPointer()

        if not QtCore.QAbstractItemModel.hasIndex(self, in_row, in_column, in_parent):
            return QtCore.QModelIndex()

        child = parent.child(in_row)
        if child:
            return QtCore.QAbstractItemModel.createIndex(self, in_row, in_column, in_parent)
        else:
            return QtCore.QModelIndex()

    def parent(self, in_index):
        if in_index.isValid():
            p = in_index.internalPointer().parent()
            if p:
                return QtCore.QAbstractItemModel.createIndex(self, p.row(),0,p)
        return QtCore.QModelIndex()

    def rowCount(self, in_index):
        if in_index.isValid():
            return in_index.internalPointer().childCount()
        return self._root.childCount()

    def columnCount(self, in_index):
        if in_index.isValid():
            return in_index.internalPointer().columnCount()
        return self._root.columnCount()

    def data(self, in_index, role):
        if not in_index.isValid():
            return None
        node = in_index.internalPointer()
        if role == QtCore.Qt.DisplayRole:
            return node.data(in_index.column())
        return None

class CustomNode(object):
    def __init__(self, in_data):
        self._data = in_data
        if type(in_data) == tuple:
            self._data = list(in_data)
        if type(in_data) in (str) or not hasattr(in_data, "__getitem__"):
            self._data = [in_data]

        self._columncount = len(self._data)
        self._children = []
        self._parent = None
        self._row = 0

    def data(self, in_column):
        if 0 <= in_column < len(self._data):
            return self._data[in_column]

    def columnCount(self):
        return self._columncount

    def childCount(self):
        return len(self._children)

    def child(self, in_row):
        if 0 <= in_row < self.childCount():
            return self._children[in_row]

    def parent(self):
        return self._parent

    def row(self):
        return self._row

    def addChild(self, in_child):
        in_child._parent = self
        in_child._row = len(self._children)
        self._children.append(in_child)
        self._columncount = max(in_child.columnCount(), self._columncount)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = Ui()
    app.exec()
