
#include "TestingDemo.h"
#include "ui_TestingDemo.h"

#include "pqEventObserver.h"
#include "pqEventSource.h"
#include "pqTestUtility.h"

#include <QApplication>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QPushButton>
#include <QTextStream>
#include <QXmlStreamAttributes>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <QtDebug>

class XMLEventObserver : public pqEventObserver
{
  QXmlStreamWriter* XMLStream;
  QString XMLString;

public:
  XMLEventObserver(QObject* p)
    : pqEventObserver(p)
  {
    this->XMLStream = NULL;
  }

  ~XMLEventObserver() override { delete this->XMLStream; }

protected:
  void setStream(QTextStream* stream) override
  {
    if (this->XMLStream)
    {
      this->XMLStream->writeEndElement();
      this->XMLStream->writeEndDocument();
      delete this->XMLStream;
      this->XMLStream = NULL;
    }
    if (this->Stream)
    {
      *this->Stream << this->XMLString;
    }
    this->XMLString = QString();
    pqEventObserver::setStream(stream);
    if (this->Stream)
    {
      this->XMLStream = new QXmlStreamWriter(&this->XMLString);
      this->XMLStream->setAutoFormatting(true);
      this->XMLStream->writeStartDocument();
      this->XMLStream->writeStartElement("events");
    }
  }

  void onRecordEvent(const QString& widget, const QString& command, const QString& arguments,
    const int& eventType) override
  {
    if (this->XMLStream)
    {
      this->XMLStream->writeStartElement("event");
      this->XMLStream->writeAttribute("widget", widget);
      if (eventType == pqEventTypes::ACTION_EVENT)
      {
        this->XMLStream->writeAttribute("command", command);
      }
      else // if(eventType == pqEventTypes::CHECK_EVENT)
      {
        this->XMLStream->writeAttribute("property", command);
      }
      this->XMLStream->writeAttribute("arguments", arguments);
      this->XMLStream->writeEndElement();
    }
  }
};

class XMLEventSource : public pqEventSource
{
  typedef pqEventSource Superclass;
  QXmlStreamReader* XMLStream;

public:
  XMLEventSource(QObject* p)
    : Superclass(p)
  {
    this->XMLStream = NULL;
  }
  ~XMLEventSource() override { delete this->XMLStream; }

protected:
  void setContent(const QString& xmlfilename) override
  {
    delete this->XMLStream;
    this->XMLStream = NULL;

    QFile xml(xmlfilename);
    if (!xml.open(QIODevice::ReadOnly))
    {
      qDebug() << "Failed to load " << xmlfilename;
      return;
    }
    QByteArray data = xml.readAll();
    this->XMLStream = new QXmlStreamReader(data);
    /* This checked for valid event objects, but also caused the first event
     * to get dropped. Commenting this out in the example. If you wish to report
     * empty XML test files a flag indicating whether valid events were found is
     * probably the best way to go.
    while (!this->XMLStream->atEnd())
      {
      QXmlStreamReader::TokenType token = this->XMLStream->readNext();
      if (token == QXmlStreamReader::StartElement)
        {
        if (this->XMLStream->name() == "event")
          {
          break;
          }
        }
      } */
    if (this->XMLStream->atEnd())
    {
      qDebug() << "Invalid xml" << endl;
    }
    return;
  }

  int getNextEvent(QString& widget, QString& command, QString& arguments, int& eventType) override
  {
    if (this->XMLStream->atEnd())
    {
      return DONE;
    }
    while (!this->XMLStream->atEnd())
    {
      QXmlStreamReader::TokenType token = this->XMLStream->readNext();
      if (token == QXmlStreamReader::StartElement)
      {
        if (this->XMLStream->name() == "event")
        {
          break;
        }
      }
    }
    if (this->XMLStream->atEnd())
    {
      return DONE;
    }
    eventType = pqEventTypes::ACTION_EVENT;
    widget = this->XMLStream->attributes().value("widget").toString();
    command = this->XMLStream->attributes().value("command").toString();
    arguments = this->XMLStream->attributes().value("arguments").toString();
    return SUCCESS;
  }
};

//-----------------------------------------------------------------------------
TestingDemo::TestingDemo()
{
  Ui::TestingDemo ui;
  ui.setupUi(this);

  QObject::connect(ui.actionRecord, SIGNAL(triggered()), this, SLOT(record()));
  QObject::connect(ui.actionPlay, SIGNAL(triggered()), this, SLOT(play()));
  QObject::connect(ui.actionOpen_Dialog, SIGNAL(triggered()), this, SLOT(popup()));

  this->TestUtility = new pqTestUtility(this);
  this->TestUtility->addEventObserver("xml", new XMLEventObserver(this));
  this->TestUtility->addEventSource("xml", new XMLEventSource(this));
}

//-----------------------------------------------------------------------------
TestingDemo::~TestingDemo()
{
  delete this->TestUtility;
}

//-----------------------------------------------------------------------------
void TestingDemo::record()
{
  QString filename =
    QFileDialog::getSaveFileName(this, "Test File Name", QString(), "XML Files (*.xml)");
  if (!filename.isEmpty())
  {
    QApplication::setActiveWindow(this);
    this->TestUtility->recordTests(filename);
  }
}

//-----------------------------------------------------------------------------
void TestingDemo::play()
{
  QString filename =
    QFileDialog::getOpenFileName(this, "Test File Name", QString(), "XML Files (*.xml)");
  if (!filename.isEmpty())
  {
    this->TestUtility->playTests(filename);
  }
}

//-----------------------------------------------------------------------------
void TestingDemo::popup()
{
  QDialog dialog;
  QHBoxLayout* hbox = new QHBoxLayout(&dialog);
  QPushButton button("Click to Close", &dialog);
  hbox->addWidget(&button);
  QObject::connect(&button, SIGNAL(clicked()), &dialog, SLOT(accept()));
  dialog.exec();
}

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  QApplication app(argc, argv);
  TestingDemo mainWindow;
  mainWindow.show();
  return app.exec();
}
