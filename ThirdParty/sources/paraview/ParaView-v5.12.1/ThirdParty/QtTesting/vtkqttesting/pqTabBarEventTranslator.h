/*=========================================================================

   Program: ParaView
   Module:    pqTabBarEventTranslator.h

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

#ifndef _pqTabBarEventTranslator_h
#define _pqTabBarEventTranslator_h

#include "pqWidgetEventTranslator.h"
#include <QPointer>
class QTabBar;

/**
\brief Records selection of a tab in a QTabBar

\note A QMainMainWindow has hidden QTabBar's that it uses to manage
QDockWidgets and other widgets that can be tabified.  These hidden QTabBars
are not named, can be numerous, and are instantiated in an unspecified order.
This makes finding the correct QTabBar difficult and sometimes impossible.
Therefore when recording selection of one of these hidden QTabBar's, the
QMainWindow is recorded as the object instead of the QTabBar itself. During
playback, if a QMainWindow was recorded then all hidden QTabBar direct
children will be searched for the tab with the recorded name. The risk with
this approach is that when multiple DockWidgets have the same windowTitle()
(which determines the tab name) then it is possible to select the wrong one.
This can be mediated in client code by ensuring that each DockWidget has a
unique windowTitle()

\sa pqEventTranslator
*/

class QTTESTING_EXPORT pqTabBarEventTranslator : public pqWidgetEventTranslator
{
  Q_OBJECT
  typedef pqWidgetEventTranslator Superclass;

public:
  pqTabBarEventTranslator(QObject* p = 0);

  using Superclass::translateEvent;
  bool translateEvent(QObject* Object, QEvent* Event, bool& Error) override;

protected Q_SLOTS:
  void indexChanged(int);

private:
  pqTabBarEventTranslator(const pqTabBarEventTranslator&);
  pqTabBarEventTranslator& operator=(const pqTabBarEventTranslator&);

  QPointer<QTabBar> CurrentObject;
};

#endif // !_pqTabBarEventTranslator_h
