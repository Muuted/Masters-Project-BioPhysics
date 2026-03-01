#from PySide6 import QtCore, QtWidgets, QtGui
#from PyQt6 import  QtCore
from PyQt6.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QVBoxLayout, QWidget, QFormLayout, QLineEdit
import sys


class AnotherWindow(QWidget):
    """
    This "window" is a QWidget. If it has no parent, it
    will appear as a free-floating window as we want.
    """
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        self.label = QLabel("Another Window")
        self.setGeometry(50,50,800,600)
        layout.addWidget(self.label)
        self.setLayout(layout)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.button1 = QPushButton("1") 
        self.input1 = QLineEdit()
        self.button2 = QPushButton("1") 
        self.input2 = QLineEdit()
        #self.setCentralWidget(self.button)
        #self.setWindowTitle("My first app")
        #self.setGeometry(50,50,800,600)

        self.layout = QFormLayout(self)
        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.input1)
        self.layout.addWidget(self.button2)
        self.layout.addWidget(self.input2)

        #self.button.clicked.connect(self.show_new_window)
        #self.setLayout(self.layout)

    def show_new_window(self,checked):
        pass#w = AnotherWindow()
        #w.show()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()