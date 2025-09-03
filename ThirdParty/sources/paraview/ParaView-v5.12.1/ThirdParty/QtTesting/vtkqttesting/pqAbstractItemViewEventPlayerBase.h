/*=========================================================================

   Program: ParaView
   Module:    pqAbstractItemViewEventPlayerBase.h

   Copyright (c) 2005,2006 Sandia Corporation, Kitware Inc.
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

========================================================================*/
#ifndef __pqAbstractItemViewEventPlayerBase_h
#define __pqAbstractItemViewEventPlayerBase_h

#include "pqWidgetEventPlayer.h"
#include <QModelIndex>

class QAbstractItemView;

/// pqAbstractItemViewEventPlayerBase is a player for QAbstractItemViewWidget. Plays back the state
/// recorded using pQAbstractItemViewEventTranslatorBase.
class QTTESTING_EXPORT pqAbstractItemViewEventPlayerBase : public pqWidgetEventPlayer
{
  Q_OBJECT
  typedef pqWidgetEventPlayer Superclass;

public:
  pqAbstractItemViewEventPlayerBase(QObject* parent = 0);
  ~pqAbstractItemViewEventPlayerBase() override = 0;

  /// Play an event on a QAbstractViewItem
  using Superclass::playEvent;
  bool playEvent(QObject* object, const QString& command, const QString& arguments, int eventType,
    bool& error) override;

protected:
  /// Get index using a previously recorded string
  static QModelIndex GetIndex(
    const QString& str_index, QAbstractItemView* abstractItemView, bool& error);

  /// Get data string using a previously recorded string
  static QString GetDataString(const QString& str_index, bool& error);

private:
  pqAbstractItemViewEventPlayerBase(const pqAbstractItemViewEventPlayerBase&); // Not implemented.
  void operator=(const pqAbstractItemViewEventPlayerBase&);                    // Not implemented.
};

#endif
