#ifndef MY_TEST_DIALOG_H
#define MY_TEST_DIALOG_H
#include "ui_mainwin.h"
class TestDialog : public QDialog
{
    Q_OBJECT
public:
    TestDialog ();
public slots:
public signals:
private:
    Ui::TestDialog ui;
};
#endif
