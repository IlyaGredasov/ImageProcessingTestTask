#include <MainWindow.hpp>
#include <QFileDialog>
#include <QMessageBox>
#include <matrix.hpp>
#include <processing.hpp>
#include <ui_MainWindow.h>

void MainWindow::openImage() {
    const QString filePath =
        QFileDialog::getOpenFileName(this, "Specify Image", "", "Images (*.png *.jpg *.jpeg *.bmp)");
    if (filePath.isEmpty()) {
        qWarning("No Image Selected");
        return;
    }
    const QPixmap pixmap(filePath);
    if (pixmap.isNull()) {
        QMessageBox::warning(this, "Error", "Failed to load image.");
        return;
    }
    ui->imageLabel->setPixmap(pixmap);
    ui->imageLabel->repaint();
}

void MainWindow::saveImage() {
    const QPixmap pixmap = ui->imageLabel->pixmap(Qt::ReturnByValue);
    if (pixmap.isNull()) {
        QMessageBox::warning(this, "Error", "No image to save.");
        return;
    }

    const QString filePath = QFileDialog::getSaveFileName(this, "Save Image", "",
        "PNG Image (*.png);;JPEG Image (*.jpg *.jpeg);;Bitmap (*.bmp)");
    if (filePath.isEmpty()) {
        return;
    }

    if (!pixmap.save(filePath)) {
        QMessageBox::critical(this, "Error", "Failed to save the image.");
    }
}

iptt::Matrix<iptt::RGB> MainWindow::asColorMatrix() const {
    const QPixmap pixmap = ui->imageLabel->pixmap(Qt::ReturnByValue);
    if (pixmap.isNull()) {
        qWarning("Pixmap is null");
        throw std::invalid_argument("Pixmap is null");
    }

    QImage image = pixmap.toImage().convertToFormat(QImage::Format_RGB888);
    const int height = image.height();
    const int width = image.width();

    iptt::Matrix<iptt::RGB> result(height, width);

    for (int i = 0; i < height; ++i) {
        const uchar* line = image.scanLine(i);
        for (int j = 0; j < width; ++j) {
            const int idx = j * 3; // RGB888 — 3 байта на пиксель
            iptt::RGB color;
            color.r = line[idx];
            color.g = line[idx + 1];
            color.b = line[idx + 2];
            result.at(i, j) = color;
        }
    }
    return result;
}

void MainWindow::fillWithColorMatrix(const iptt::Matrix<iptt::RGB>& m) const {
    const auto height = static_cast<int>(m.rows());
    const auto width = static_cast<int>(m.cols());
    QImage image(width, height, QImage::Format_RGB888);

    for (int i = 0; i < height; ++i) {
        uchar* line = image.scanLine(i);
        for (int j = 0; j < width; ++j) {
            const auto& [r, g, b] = m.at(i, j);
            const int idx = j * 3;
            line[idx + 0] = r;
            line[idx + 1] = g;
            line[idx + 2] = b;
        }
    }
    ui->imageLabel->setPixmap(QPixmap::fromImage(image));
    ui->imageLabel->repaint();
}

iptt::Matrix<int> MainWindow::asGrayscaleMatrix() const {
    const QPixmap pixmap = ui->imageLabel->pixmap(Qt::ReturnByValue);
    if (pixmap.isNull()) {
        qWarning("Pixmap is null");
        throw std::invalid_argument("Pixmap is null");
    }

    QImage image = pixmap.toImage().convertToFormat(QImage::Format_RGB888);
    const int height = image.height();
    const int width = image.width();

    iptt::Matrix<int> result(height, width);

    for (int i = 0; i < height; ++i) {
        const uchar* line = image.scanLine(i);
        for (int j = 0; j < width; ++j) {
            const int idx = j * 3;
            const auto gray = iptt::brightness(iptt::RGB{line[idx + 0], line[idx + 1], line[idx + 2]});
            result.at(i, j) = static_cast<int>(gray);
        }
    }
    return result;
}

void MainWindow::fillWithGrayscaleMatrix(const iptt::Matrix<int>& m) const {
    const auto height = static_cast<int>(m.rows());
    const auto width = static_cast<int>(m.cols());

    QImage image(width, height, QImage::Format_Grayscale8);

    for (int i = 0; i < height; ++i) {
        uchar* line = image.scanLine(i);
        for (int j = 0; j < width; ++j) {
            line[j] = std::clamp(m.at(i, j), 0, 255);
        }
    }

    ui->imageLabel->setPixmap(QPixmap::fromImage(image));
    ui->imageLabel->repaint();
}

void MainWindow::applyGaussianFilter() const { fillWithColorMatrix(iptt::gaussian_filter(asColorMatrix())); }
void MainWindow::applySobelOperator() const { fillWithGrayscaleMatrix(iptt::sobel_operator(asGrayscaleMatrix())); }

void MainWindow::applyErosion() const {
    const auto thresh = ui->firstLineEdit->text().toInt();
    fillWithGrayscaleMatrix(iptt::erose(asGrayscaleMatrix(), thresh));
}

void MainWindow::applyDilatation() const {
    const auto thresh = ui->firstLineEdit->text().toInt();
    fillWithGrayscaleMatrix(iptt::dilate(asGrayscaleMatrix(), thresh));
}

void MainWindow::applyGrayscale() const { fillWithGrayscaleMatrix(iptt::grayscale(asColorMatrix())); }
void MainWindow::applyContrastStretch() const { fillWithColorMatrix(iptt::contrast_stretch(asColorMatrix())); }

void MainWindow::applyScale() const {
    const auto sx = ui->firstLineEdit->text().toDouble();
    const auto sy = ui->secondLineEdit->text().toDouble();
    fillWithColorMatrix(iptt::scale(asColorMatrix(), sx, sy));
}

void MainWindow::applyFlip() const {
    const auto h = ui->firstLineEdit->text().toInt();
    const auto v = ui->secondLineEdit->text().toInt();
    fillWithColorMatrix(iptt::flip(asColorMatrix(), h != 0, v != 0));
}

void MainWindow::applyRotate() const {
    const auto pi = std::atan(1) * 4;
    const auto angle = ui->firstLineEdit->text().toDouble() * pi / 180.0;
    const auto rotated = iptt::rotate(asColorMatrix(), angle);
    const auto cropped = iptt::crop_borders(rotated, [](const iptt::RGB& rgb) { return iptt::brightness(rgb) > 50; });
    fillWithColorMatrix(cropped);
}

void MainWindow::applyTransformComponent() const {
    const auto r = ui->firstLineEdit->text().toShort();
    const auto g = ui->secondLineEdit->text().toShort();
    const auto b = ui->thirdLineEdit->text().toShort();
    fillWithColorMatrix(iptt::transform(asColorMatrix(), [r, b, g](const iptt::RGB& cell) {
        return iptt::RGB{iptt::clamp(cell.r + r), iptt::clamp(cell.g + g), iptt::clamp(cell.b + b)};
    }));
}

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent), ui(new Ui_MainWindow) {
    ui->setupUi(this);
    connect(ui->openFileAction, &QAction::triggered, this, &MainWindow::openImage);
    connect(ui->saveFileAction, &QAction::triggered, this, &MainWindow::saveImage);
    connect(ui->gaussButton, &QPushButton::clicked, this, &MainWindow::applyGaussianFilter);
    connect(ui->sobelButton, &QPushButton::clicked, this, &MainWindow::applySobelOperator);
    connect(ui->erosionButton, &QPushButton::clicked, this, &MainWindow::applyErosion);
    connect(ui->dilatationButton, &QPushButton::clicked, this, &MainWindow::applyDilatation);
    connect(ui->grayscaleButton, &QPushButton::clicked, this, &MainWindow::applyGrayscale);
    connect(ui->contrastButton, &QPushButton::clicked, this, &MainWindow::applyContrastStretch);
    connect(ui->scaleButton, &QPushButton::clicked, this, &MainWindow::applyScale);
    connect(ui->flipButton, &QPushButton::clicked, this, &MainWindow::applyFlip);
    connect(ui->rotateButton, &QPushButton::clicked, this, &MainWindow::applyRotate);
    connect(ui->transformComponentButton, &QPushButton::clicked, this, &MainWindow::applyTransformComponent);
}

MainWindow::~MainWindow() { delete ui; }
