/*=========================================================================

   Program: ParaView
   Module:    pqWidgetEventPlayer.cxx

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

#include "pqWidgetEventPlayer.h"

#include <QContextMenuEvent>
#include <QCoreApplication>
#include <QPoint>
#include <QWidget>
#include <QtDebug>
#include <cmath> // for std::abs

#include "pqEventTypes.h"

pqWidgetEventPlayer::pqWidgetEventPlayer(QObject* p)
  : QObject(p)
{
}

pqWidgetEventPlayer::~pqWidgetEventPlayer()
{
}

bool pqWidgetEventPlayer::playEvent(
  QObject* object, const QString& command, const QString& arguments, bool& error)
{
  QWidget* widget = qobject_cast<QWidget*>(object);
  if (widget)
  {
    if (command == "contextMenu")
    {
      QPoint pt(widget->x(), widget->y());
      QPoint globalPt = widget->mapToGlobal(pt);
      QContextMenuEvent e(QContextMenuEvent::Other, pt, globalPt);
      qApp->notify(widget, &e);
      return true;
    }
    else if (command == "size")
    {
      QStringList entries = arguments.split(',');
      if (entries.size() == 2)
      {
        QSize sz = widget->size();
        error = (sz.width() != entries[0].toInt() || sz.height() != entries[1].toInt());
        if (error)
        {
          qCritical() << "Size mismatch: (" << arguments << ") != " << sz;
        }
        return true;
      }
    }
    else if (command == "sizeGreaterOrEqual")
    {
      QStringList entries = arguments.split(',');
      if (entries.size() == 2)
      {
        QSize sz = widget->size();
        error = (sz.width() < entries[0].toInt() || sz.height() < entries[1].toInt());
        if (error)
        {
          qCritical() << "Size mismatch: (" << arguments << ") > " << sz;
        }
        return true;
      }
    }
    else if (command == "aspectRatio")
    {
      double targetAspectRatio = arguments.toDouble();
      QSize sz = widget->size();
      double aspectRatio = static_cast<double>(sz.width()) / sz.height();
      error = (std::abs(targetAspectRatio - aspectRatio) > 1e-2);
      if (error)
      {
        qCritical() << "Ratio mismatch: (" << targetAspectRatio << ") != " << aspectRatio;
      }
      return true;
    }
  }
  return false;
}

bool pqWidgetEventPlayer::playEvent(
  QObject* object, const QString& command, const QString& arguments, int eventType, bool& error)
{
  if (eventType == pqEventTypes::ACTION_EVENT)
  {
    return this->playEvent(object, command, arguments, error);
  }
  return false;
}
