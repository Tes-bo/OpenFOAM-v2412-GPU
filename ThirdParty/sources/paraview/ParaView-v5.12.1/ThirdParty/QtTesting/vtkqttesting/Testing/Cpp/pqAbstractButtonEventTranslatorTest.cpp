/*=========================================================================

   Program: ParaView

   Copyright (c) 2005-2008 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2.

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

// Qt includes
#include <QAction>
#include <QApplication>
#include <QMenu>
#include <QToolButton>
#include <QtTest/QtTest>

// QtTesting includes
#include "pqTestUtility.h"

#include "pqTest.h"

// ----------------------------------------------------------------------------
class pqAbstractButtonEventTranslatorTester : public QObject
{
  Q_OBJECT

private Q_SLOTS:

  void initTestCase();
  void init();
  void cleanup();
  void cleanupTestCase();

  void testToolButton();
  void testToolButton_data();

private:
  QToolButton* ToolButton;

  pqTestUtility* TestUtility;
  pqDummyEventObserver* EventObserver;
};

// ----------------------------------------------------------------------------
void pqAbstractButtonEventTranslatorTester::initTestCase()
{
  this->TestUtility = new pqTestUtility();
  this->TestUtility->setRecordWithDialog(false);
  this->EventObserver = new pqDummyEventObserver();
  this->TestUtility->addEventObserver("xml", this->EventObserver);
  this->TestUtility->addEventSource("xml", new pqDummyEventSource());

  this->ToolButton = 0;
}

// ----------------------------------------------------------------------------
void pqAbstractButtonEventTranslatorTester::init()
{
  // Init the DoubleSpinBox
  this->ToolButton = new QToolButton();
  this->ToolButton->setText("toolButton");
  this->ToolButton->setObjectName("toolButton");
  this->ToolButton->resize(200, 200);

  this->ToolButton->show();
#if QT_VERSION < QT_VERSION_CHECK(5, 3, 0)
  QTest::qWaitForWindowShown(this->ToolButton);
#else
  QTest::qWaitForWindowExposed(this->ToolButton);
#endif

  // Start to record events
  this->TestUtility->recordTestsBySuffix("xml");

  // Fire the event "enter" to connect ToolButton signals to the translator slots
  QEvent enter(QEvent::Enter);
  qApp->notify(this->ToolButton, &enter);
}

// ----------------------------------------------------------------------------
void pqAbstractButtonEventTranslatorTester::cleanup()
{
  this->EventObserver->Text = QString();
  this->TestUtility->stopRecords(0);
  delete this->ToolButton;
  this->ToolButton = 0;
}

// ----------------------------------------------------------------------------
void pqAbstractButtonEventTranslatorTester::cleanupTestCase()
{
  this->EventObserver = 0;

  delete this->TestUtility;
}

// ----------------------------------------------------------------------------
void pqAbstractButtonEventTranslatorTester::testToolButton()
{
  QFETCH(int, popupMode);
  this->ToolButton->setPopupMode(static_cast<QToolButton::ToolButtonPopupMode>(popupMode));

  QFETCH(bool, withDefaultAction);
  QAction* defaultAction = 0;
  if (withDefaultAction)
  {
    defaultAction = new QAction("action", this->ToolButton);
    defaultAction->setObjectName("action");
    this->ToolButton->setDefaultAction(defaultAction);
  }

  QFETCH(bool, checkable);
  if (defaultAction)
  {
    defaultAction->setCheckable(checkable);
  }
  else
  {
    this->ToolButton->setCheckable(checkable);
  }

  QFETCH(bool, withMenu);
  if (withMenu)
  {
    QMenu* menu = new QMenu(this->ToolButton);
    QAction* menuAction = new QAction("menu action", menu);
    menu->addAction(menuAction);
    if (defaultAction)
    {
      defaultAction->setMenu(menu);
    }
    else
    {
      this->ToolButton->setMenu(menu);
    }
    // After opening the menu, the event loop might get blocked.
    QTimer::singleShot(1000, menu, SLOT(close()));
  }

  QTest::mousePress(this->ToolButton, Qt::LeftButton, Qt::NoModifier);
  QFETCH(bool, longClick);
  QTest::mouseRelease(
    this->ToolButton, Qt::LeftButton, Qt::NoModifier, QPoint(), longClick ? 800 : -1);

  QFETCH(QString, recordEmitted);
  bool res = this->EventObserver->Text == recordEmitted;
  QCOMPARE(this->EventObserver->Text, recordEmitted);
}

// ----------------------------------------------------------------------------
void pqAbstractButtonEventTranslatorTester::testToolButton_data()
{
  QTest::addColumn<int>("popupMode");
  QTest::addColumn<bool>("withDefaultAction");
  QTest::addColumn<bool>("checkable");
  QTest::addColumn<bool>("withMenu");
  QTest::addColumn<bool>("longClick");
  QTest::addColumn<QString>("recordEmitted");

  // QTest::newRow("fail") << 2 << false << false << true << false<< "toolButton, activate, #";
  for (int i = 0; i < 0x20; ++i)
  {
    int popup = i & 0x01 ? QToolButton::DelayedPopup : QToolButton::InstantPopup;
    const bool withDefaultAction = i & 0x02;
    const bool checkable = i & 0x04;
    const bool withMenu = i & 0x08;
    const bool longClick = i & 0x10;

    QString testName =
      QString("popup=%1 withDefaultAction=%2 checkable=%3 withMenu=%4 longClick=%5")
        .arg(popup)
        .arg(withDefaultAction)
        .arg(checkable)
        .arg(withMenu)
        .arg(longClick);

    bool longActivate = withMenu && popup == QToolButton::DelayedPopup && longClick;
    QString recordEmitted =
      QString("toolButton%1, %2, %3#")
        .arg(withDefaultAction && !longActivate ? "/action" : "")
        .arg(longActivate ? "longActivate" : (checkable ? "set_boolean" : "activate"))
        .arg(checkable && !longActivate ? "true" : "");

    QTest::newRow(testName.toUtf8()) << popup << withDefaultAction << checkable << withMenu
                                     << longClick << recordEmitted;
  }
}

// ----------------------------------------------------------------------------
CTK_TEST_MAIN(pqAbstractButtonEventTranslatorTest)
#include "moc_pqAbstractButtonEventTranslatorTest.cpp"
