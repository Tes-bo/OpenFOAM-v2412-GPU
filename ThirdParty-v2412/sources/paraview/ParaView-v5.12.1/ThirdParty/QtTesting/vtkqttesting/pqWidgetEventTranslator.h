/*=========================================================================

   Program: ParaView
   Module:    pqWidgetEventTranslator.h

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

#ifndef _pqWidgetEventTranslator_h
#define _pqWidgetEventTranslator_h

#include "QtTestingExport.h"
#include <QObject>
#include <QRect>

/**
Abstract interface for an object that can translate
low-level Qt events into high-level, serializable ParaView
events, for test-cases, demos, tutorials, etc.

\sa pqEventTranslator
*/
class QTTESTING_EXPORT pqWidgetEventTranslator : public QObject
{
  Q_OBJECT

  typedef QObject Superclass;

public:
  pqWidgetEventTranslator(QObject* p = 0);
  ~pqWidgetEventTranslator() override;

  /** Derivatives should implement this and translate events into commands,
  returning "true" if they handled the event, and setting Error
  to "true" if there were any problems. eventType allow to specify different types of events
  like check event*/
  virtual bool translateEvent(QObject* object, QEvent* event, bool& error);
  virtual bool translateEvent(QObject* object, QEvent* event, int eventType, bool& error);

Q_SIGNALS:
  /// Derivatives should emit this signal whenever they wish to record a high-level event
  void recordEvent(
    QObject* Object, const QString& Command, const QString& Arguments, int eventType);
  void recordEvent(QObject* Object, const QString& Command, const QString& Arguments);
  void specificOverlay(const QRect& geometry);

protected:
  pqWidgetEventTranslator(const pqWidgetEventTranslator&);
  pqWidgetEventTranslator& operator=(const pqWidgetEventTranslator&);
};

#endif // !_pqWidgetEventTranslator_h
