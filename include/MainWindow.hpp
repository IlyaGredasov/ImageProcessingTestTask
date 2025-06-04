#pragma once

#include <QMainWindow>
#include <matrix.hpp>
#include <processing.hpp>

class Ui_MainWindow;

class MainWindow final : public QMainWindow {
	Q_OBJECT

public:
	explicit MainWindow(QWidget* parent = nullptr);
	iptt::Matrix<iptt::RGB> asColorMatrix() const;
	void fillWithColorMatrix(const iptt::Matrix<iptt::RGB>& m) const;
	iptt::Matrix<int> asGrayscaleMatrix() const;
	void fillWithGrayscaleMatrix(const iptt::Matrix<int>& m) const;

	void applyGaussianFilter() const;
	void applySobelOperator() const;
	void applyErosion() const;
	void applyDilatation() const;
	void applyGrayscale() const;
	void applyContrastStretch() const;
	void applyScale() const;
	void applyFlip() const;
	void applyRotate() const;
	void applyTransformComponent() const;
	~MainWindow() override;

private slots:
	void openImage();
	void saveImage();

private:
	Ui_MainWindow* ui;
};
