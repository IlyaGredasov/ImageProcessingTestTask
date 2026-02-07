#include "MainWindow.hpp"
#include <QApplication>

int main(int argc, char* argv[]) {
    const QApplication app(argc, argv);
    MainWindow window;
    window.adjustSize();
    window.show();
    return QApplication::exec();
}
