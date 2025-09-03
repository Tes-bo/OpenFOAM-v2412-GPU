/*=========================================================================

   Program: ParaView
   Module:    pqBasicWidgetEventPlayer.cxx

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

#include "pqBasicWidgetEventPlayer.h"

#include <QApplication>
#include <QKeyEvent>
#include <QWidget>
#include <QtDebug>

#include "pqEventTypes.h"

pqBasicWidgetEventPlayer::pqBasicWidgetEventPlayer(QObject* p)
  : pqWidgetEventPlayer(p)
{
}

bool pqBasicWidgetEventPlayer::playEvent(
  QObject* object, const QString& command, const QString& arguments, int eventType, bool& error)
{
  QWidget* widget = qobject_cast<QWidget*>(object);
  if (widget)
  {
    if (eventType == pqEventTypes::ACTION_EVENT)
    {
      {
        if (command == "key")
        {
          QKeyEvent kd(QEvent::KeyPress, arguments.toInt(), Qt::NoModifier);
          QKeyEvent ku(QEvent::KeyRelease, arguments.toInt(), Qt::NoModifier);
          qApp->notify(widget, &kd);
          qApp->notify(widget, &ku);
          return true;
        }
        else if (command == "keyEvent")
        {
          // The double split is for retrocompatibility
          QStringList data = arguments.split(',');
          if (data.size() != 6)
          {
            data = arguments.split(':');
          }

          if (data.size() == 6)
          {
            QKeyEvent ke(static_cast<QEvent::Type>(data[0].toInt()), data[1].toInt(),
              static_cast<Qt::KeyboardModifiers>(data[2].toInt()), data[3], !!data[4].toInt(),
              data[5].toInt());
            qApp->notify(widget, &ke);
            return true;
          }
        }
        else if (command.startsWith("mouse"))
        {
          QStringList args = arguments.split(',');
          if (args.size() == 5)
          {
            Qt::MouseButtons buttons = static_cast<Qt::MouseButton>(args[1].toInt());
            Qt::KeyboardModifiers keym = static_cast<Qt::KeyboardModifier>(args[2].toInt());
            int x = args[3].toInt();
            int y = args[4].toInt();
            QPoint pt(x, y);
            if (command == "mouseWheel")
            {
              int delta = args[0].toInt();
              QWheelEvent we(QPoint(x, y), QPoint(x, y), QPoint(0, 0), QPoint(0, delta), buttons,
                keym, Qt::NoScrollPhase, false);
              QCoreApplication::sendEvent(object, &we);
              return true;
            }
            Qt::MouseButton button = static_cast<Qt::MouseButton>(args[0].toInt());
            QEvent::Type type = QEvent::MouseButtonPress;
            type = command == "mouseMove" ? QEvent::MouseMove : type;
            type = command == "mouseRelease" ? QEvent::MouseButtonRelease : type;
            type = command == "mouseDblClick" ? QEvent::MouseButtonDblClick : type;
            if (type == QEvent::MouseMove)
            {
              // We have not been setting mouse move correctly.
              buttons = button;
              button = Qt::NoButton;
            }
            QMouseEvent e(type, pt, button, buttons, keym);
            qApp->notify(widget, &e);
            return true;
          }
        }
      }
    }
    else if (eventType == pqEventTypes::CHECK_EVENT)
    {
      // Recover QProperty
      QVariant propertyValue = object->property(command.toUtf8().data());

      // Check it is valid
      if (!propertyValue.isValid())
      {
        QString errorMessage = object->objectName() + " has no valid property named:" + command;
        qCritical() << errorMessage.toUtf8().data();
        error = true;
        return true;
      }

      // Check property value
      if (propertyValue.toString().replace("\t", " ") != arguments)
      {
        QString errorMessage = object->objectName() + " property value is: " +
          propertyValue.toString() + ". Expecting: " + arguments + ".";
        qCritical() << errorMessage.toUtf8().data();
        error = true;
      }
      return true;
    }
  }
  return this->Superclass::playEvent(object, command, arguments, error);
}
