from PyQt5 import QtWidgets, uic, QtCore
from atoms import Atom, EnergyLevel, Transition
from atom_import import pickle_atom, load_from_pickle
from atom_library import populate_levels
from parsers import parse_NIST_levels
import sys
import numpy as np
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
    def __init__(self, nodes):
        super(AtomModel, self).__init__()
        self._root = CustomNode(None)
        print(nodes)
        for node in nodes:
            self._root.addChild(node)

    def addChild(self, node, _parent):
        if not _parent or not _parent.isValid():
            parent = self._root
        else:
            parent = _parent.internalPointer()
        parent.addChild(node)

    def index(self, row, column, _parent=None):
        if not _parent or not _parent.isValid():
            parent = self._root
        else:
            parent = _parent.internalPointer()

        if not QtCore.QAbstractItemModel.hasIndex(self, row, column, _parent):
            return QtCore.QModelIndex()

        child = parent.child(row)
        if child:
            return QtCore.QAbstractItemModel.createIndex(self, row, column, child)
        else:
            return QtCore.QModelIndex()

    def parent(self, index):
        if index.isValid():
            p = index.internalPointer().parent()
            if p:
                return QtCore.QAbstractItemModel.createIndex(self, p.row(), 0, p)
        return QtCore.QModelIndex()

    def rowCount(self, index):
        if index.isValid():
            return index.internalPointer().childCount()
        return self._root.childCount()

    def columnCount(self, index):
        if index.isValid():
            return index.internalPointer().columnCount()
        return self._root.columnCount()

    def data(self, index, role):
        if not index.isValid():
            return None
        node = index.internalPointer()
        if role == QtCore.Qt.DisplayRole:
            return node.data(index.column())
        return None

class CustomNode(object):
    def __init__(self, data):
        self._data = data
        if type(data) == tuple:
            self._data = list(data)
        if type(data) == EnergyLevel:
            self._data = [data]
        if type(data) is str or not hasattr(data, "__getitem__"):
            self._data = [data]

        self._columncount = len(self._data)
        self._children = []
        self._parent = None
        self._row = 0

    def data(self, column):
        if 0 <= column < len(self._data):
            if type(self._data[column]) == EnergyLevel:
                return self._data[column].name
            else:
                return self._data[column]

    def columnCount(self):
        return self._columncount

    def childCount(self):
        return len(self._children)

    def child(self, row):
        if 0 <= row < self.childCount():
            return self._children[row]

    def parent(self):
        return self._parent

    def row(self):
        return self._row

    def addChild(self, child):
        child._parent = self
        child._row = len(self._children)
        self._children.append(child)
        self._columncount = max(child.columnCount(), self._columncount)

class AtomTree:
    def __init__(self, atom):
        self.items = []
        self.atom = atom
        for name, level in self.atom.levels.items():
            print(name)
            print(level)
            print(level.Fs)
            l = CustomNode(level)
            for F in level.Fs:
                f = CustomNode([F])
                for mf in np.arange(-F, F+1, 1):
                    f.addChild(CustomNode([str(mf)]))
                l.addChild(f)
            self.items.append(l)

        self.tw = QtWidgets.QTreeView()
        self.tw.setModel(AtomModel(self.items))


if __name__ == "__main__":
    Yb171 = load_from_pickle("C:/Users/jippi/PycharmProjects/Atom_Lines_Tool/171Yb.atom")

    app = QtWidgets.QApplication(sys.argv)
    # window = Ui()
    # app.exec()

    myAtom = AtomTree(Yb171)
    myAtom.tw.show()
    app.exec()
