# =============================
# File: mtxviewer/style.py
# =============================
STYLE = r"""
/* ---------- Base ---------- */
* { outline: 0; }
QWidget {
    background-color: #2b2d31;
    color: #dcddde;
    font-family: 'Inter', 'Helvetica Neue', Arial, sans-serif;
    font-size: 15px;
}

/* Panels / group boxes */
QToolBar { background: #1e1f22; border-bottom: 1px solid #202225; }
QGroupBox {
    border: 1px solid #3a3c43;
    border-radius: 8px;
    margin-top: 10px;
    background: #313338;
}
QGroupBox::title {
    left: 10px;
    padding: 0 6px;
    color: #b9bbbe;
    background: transparent;
}
QSplitter::handle { background: #202225; }

/* ---------- Buttons ---------- */
QPushButton[primary="true"] {
    background-color: #5865f2;  /* Discord blurple */
    color: #ffffff;
    border: none;
    border-radius: 6px;
    padding: 8px 14px;
    font-weight: 600;
}
QPushButton[primary="true"]:hover { background-color: #4752c4; }
QPushButton[primary="true"]:disabled { background-color: #40444b; color: #72767d; }

QPushButton {
    background: #3a3c43;
    border: 1px solid #202225;
    border-radius: 6px;
    padding: 6px 12px;
    color: #dcddde;
}
QPushButton:hover { background: #4e5058; }

QToolButton {
    background: transparent;
    border: none;
    padding: 6px;
    border-radius: 6px;
    color: #b9bbbe;
}
QToolButton:hover {
    background: rgba(255,255,255,0.07);
    color: #ffffff;
}

/* ---------- Inputs ---------- */
QComboBox, QSpinBox, QLineEdit {
    background-color: #1e1f22;
    border: 1px solid #3a3c43;
    border-radius: 6px;
    padding: 6px 10px;
    min-height: 28px;
    color: #dcddde;
}
QComboBox QAbstractItemView {
    background: #2b2d31;
    color: #dcddde;
    border: 1px solid #3a3c43;
    selection-background-color: #5865f2;
    selection-color: #ffffff;
}
QLabel { background: transparent; }

/* ---------- Tables / Views ---------- */
QTableWidget, QTreeView, QListView {
    background: #1e1f22;
    border: 1px solid #3a3c43;
    border-radius: 6px;
    gridline-color: #3a3c43;
    selection-background-color: #4752c4;
    selection-color: #ffffff;
}
QHeaderView::section {
    background: #2b2d31;
    color: #dcddde;
    border: none;
    padding: 8px 10px;
    font-weight: 600;
    font-size: 14px;
}
QTableCornerButton::section {
    background: #2b2d31;
    border: none;
}

/* ---------- Tabs ---------- */
QTabWidget::pane {
    border: 1px solid #3a3c43;
    border-radius: 6px;
    background: #313338;
}
QTabBar::tab {
    background: #2b2d31;
    color: #b9bbbe;
    padding: 8px 14px;
    margin-right: 4px;
    border: 1px solid #3a3c43;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    min-width: 120px;
    font-weight: 600;
}
QTabBar::tab:hover {
    background: #3a3c43;
    color: #ffffff;
}
QTabBar::tab:selected {
    background: #313338;
    color: #ffffff;
    border-bottom: 2px solid #5865f2;
}

/* ---------- Scrollbars ---------- */
QScrollBar:vertical, QScrollBar:horizontal {
    background: #2b2d31;
    border: none;
    margin: 0;
}
QScrollBar:vertical { width: 10px; }
QScrollBar:horizontal { height: 10px; }
QScrollBar::handle:vertical, QScrollBar::handle:horizontal {
    background: #5865f2;
    border-radius: 5px;
    min-height: 20px; min-width: 20px;
}
QScrollBar::handle:vertical:hover, QScrollBar::handle:horizontal:hover {
    background: #4752c4;
}
QScrollBar::add-line, QScrollBar::sub-line,
QScrollBar::add-page, QScrollBar::sub-page {
    background: transparent; border: none;
}

/* ---------- Text edits ---------- */
QPlainTextEdit, QTextEdit {
    background: #1e1f22;
    border: 1px solid #3a3c43;
    border-radius: 6px;
    color: #dcddde;
}
"""
