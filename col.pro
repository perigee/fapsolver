# -------------------------------------------------
# Project created by QtCreator 2009-10-05T14:27:11
# -------------------------------------------------
QT -= core \
    gui

# CONFIG += link_pkgconfig
# PKGCONFIG += libxml-2.0
TARGET = solver
CONFIG += console
CONFIG -= app_bundle
TEMPLATE = app
SOURCES += main.cpp
HEADERS += dataIO.hpp \
    colAlgo.hpp \
    Base.hpp \
    wave.hpp \
    heuristic.hpp \
    statistics.hpp \
    algorithm.hpp \
    heuristic.hpp \
    color.hpp
