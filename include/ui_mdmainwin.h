/********************************************************************************
** Form generated from reading UI file 'mdmainwin.ui'
**
** Created: Mon 25. Mar 23:27:18 2013
**      by: Qt User Interface Compiler version 4.8.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MDMAINWIN_H
#define UI_MDMAINWIN_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTextBrowser>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_mdmainwin
{
public:
    QAction *actionExit;
    QAction *actionAbout;
    QAction *actionSave_Settings;
    QAction *actionLoad_Settings;
    QWidget *centralWidget;
    QVBoxLayout *verticalLayout_2;
    QTabWidget *tabWidget_2;
    QWidget *tab_settings;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QWidget *tab_simulation_2;
    QVBoxLayout *verticalLayout_4;
    QHBoxLayout *horizontalLayout_3;
    QPushButton *start_simulation_pb;
    QCheckBox *close_when_finished_cb;
    QSpacerItem *horizontalSpacer_2;
    QHBoxLayout *horizontalLayout_2;
    QTextBrowser *simulation_output_tb;
    QWidget *tab_plots;
    QWidget *tab_GL;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout_6;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuHelp;
    QMenu *menuSettings;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *mdmainwin)
    {
        if (mdmainwin->objectName().isEmpty())
            mdmainwin->setObjectName(QString::fromUtf8("mdmainwin"));
        mdmainwin->resize(1079, 619);
        QSizePolicy sizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(mdmainwin->sizePolicy().hasHeightForWidth());
        mdmainwin->setSizePolicy(sizePolicy);
        actionExit = new QAction(mdmainwin);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionAbout = new QAction(mdmainwin);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        actionSave_Settings = new QAction(mdmainwin);
        actionSave_Settings->setObjectName(QString::fromUtf8("actionSave_Settings"));
        actionLoad_Settings = new QAction(mdmainwin);
        actionLoad_Settings->setObjectName(QString::fromUtf8("actionLoad_Settings"));
        centralWidget = new QWidget(mdmainwin);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        verticalLayout_2 = new QVBoxLayout(centralWidget);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        tabWidget_2 = new QTabWidget(centralWidget);
        tabWidget_2->setObjectName(QString::fromUtf8("tabWidget_2"));
        tab_settings = new QWidget();
        tab_settings->setObjectName(QString::fromUtf8("tab_settings"));
        verticalLayout = new QVBoxLayout(tab_settings);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));

        verticalLayout->addLayout(horizontalLayout);

        tabWidget_2->addTab(tab_settings, QString());
        tab_simulation_2 = new QWidget();
        tab_simulation_2->setObjectName(QString::fromUtf8("tab_simulation_2"));
        verticalLayout_4 = new QVBoxLayout(tab_simulation_2);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setContentsMargins(11, 11, 11, 11);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        start_simulation_pb = new QPushButton(tab_simulation_2);
        start_simulation_pb->setObjectName(QString::fromUtf8("start_simulation_pb"));

        horizontalLayout_3->addWidget(start_simulation_pb);

        close_when_finished_cb = new QCheckBox(tab_simulation_2);
        close_when_finished_cb->setObjectName(QString::fromUtf8("close_when_finished_cb"));
        close_when_finished_cb->setEnabled(true);

        horizontalLayout_3->addWidget(close_when_finished_cb);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_2);


        verticalLayout_4->addLayout(horizontalLayout_3);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));

        verticalLayout_4->addLayout(horizontalLayout_2);

        simulation_output_tb = new QTextBrowser(tab_simulation_2);
        simulation_output_tb->setObjectName(QString::fromUtf8("simulation_output_tb"));

        verticalLayout_4->addWidget(simulation_output_tb);

        tabWidget_2->addTab(tab_simulation_2, QString());
        tab_plots = new QWidget();
        tab_plots->setObjectName(QString::fromUtf8("tab_plots"));
        tabWidget_2->addTab(tab_plots, QString());
        tab_GL = new QWidget();
        tab_GL->setObjectName(QString::fromUtf8("tab_GL"));
        sizePolicy.setHeightForWidth(tab_GL->sizePolicy().hasHeightForWidth());
        tab_GL->setSizePolicy(sizePolicy);
        verticalLayoutWidget = new QWidget(tab_GL);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(10, 10, 1031, 501));
        verticalLayout_6 = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setContentsMargins(11, 11, 11, 11);
        verticalLayout_6->setObjectName(QString::fromUtf8("verticalLayout_6"));
        verticalLayout_6->setContentsMargins(0, 0, 0, 0);
        tabWidget_2->addTab(tab_GL, QString());

        verticalLayout_2->addWidget(tabWidget_2);

        mdmainwin->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(mdmainwin);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1079, 26));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuHelp = new QMenu(menuBar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        menuSettings = new QMenu(menuBar);
        menuSettings->setObjectName(QString::fromUtf8("menuSettings"));
        mdmainwin->setMenuBar(menuBar);
        statusbar = new QStatusBar(mdmainwin);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        mdmainwin->setStatusBar(statusbar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuSettings->menuAction());
        menuBar->addAction(menuHelp->menuAction());
        menuFile->addAction(actionExit);
        menuHelp->addAction(actionAbout);
        menuSettings->addAction(actionSave_Settings);
        menuSettings->addAction(actionLoad_Settings);

        retranslateUi(mdmainwin);
        QObject::connect(actionExit, SIGNAL(triggered()), mdmainwin, SLOT(close()));

        tabWidget_2->setCurrentIndex(3);


        QMetaObject::connectSlotsByName(mdmainwin);
    } // setupUi

    void retranslateUi(QMainWindow *mdmainwin)
    {
        mdmainwin->setWindowTitle(QApplication::translate("mdmainwin", "meshworks", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("mdmainwin", "Exit", 0, QApplication::UnicodeUTF8));
        actionAbout->setText(QApplication::translate("mdmainwin", "About...", 0, QApplication::UnicodeUTF8));
        actionSave_Settings->setText(QApplication::translate("mdmainwin", "Save Settings...", 0, QApplication::UnicodeUTF8));
        actionLoad_Settings->setText(QApplication::translate("mdmainwin", "Load Settings...", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_settings), QApplication::translate("mdmainwin", "Settings", 0, QApplication::UnicodeUTF8));
        start_simulation_pb->setText(QApplication::translate("mdmainwin", "Start Simulation", 0, QApplication::UnicodeUTF8));
        close_when_finished_cb->setText(QApplication::translate("mdmainwin", "Close when finished", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_simulation_2), QApplication::translate("mdmainwin", "Simulation", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_plots), QApplication::translate("mdmainwin", "Plots", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_GL), QApplication::translate("mdmainwin", "GL", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("mdmainwin", "File", 0, QApplication::UnicodeUTF8));
        menuHelp->setTitle(QApplication::translate("mdmainwin", "Help", 0, QApplication::UnicodeUTF8));
        menuSettings->setTitle(QApplication::translate("mdmainwin", "Settings", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class mdmainwin: public Ui_mdmainwin {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MDMAINWIN_H
