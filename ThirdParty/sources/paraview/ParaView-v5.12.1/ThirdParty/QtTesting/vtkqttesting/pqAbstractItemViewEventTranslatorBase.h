/*=========================================================================

   Program: ParaView
   Module:    pqAbstractItemViewEventTranslatorBase.h

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
#ifndef __pqAbstractItemViewEventTranslatorBase_h
#define __pqAbstractItemViewEventTranslatorBase_h

#include "pqWidgetEventTranslator.h"
#include <QAbstractItemModel> // for QModelIndexList
#include <QPointer>

class QModelIndex;
class QAbstractItemView;
class QItemSelection;
class QItemSelectionModel;

/// Event recorder for QAbstractItemView. Records the toggling of the check states for
/// tree widget items. The recorded state can be played back using
/// pqTreeViewEventPlayer.
class QTTESTING_EXPORT pqAbstractItemViewEventTranslatorBase : public pqWidgetEventTranslator
{
  Q_OBJECT
  typedef pqWidgetEventTranslator Superclass;

public:
  pqAbstractItemViewEventTranslatorBase(QObject* parent = 0);
  ~pqAbstractItemViewEventTranslatorBase() override;

  /// Overridden to handle events on QAbstractItemView
  using pqWidgetEventTranslator::translateEvent;
  bool translateEvent(QObject* object, QEvent* event, int eventType, bool& error) override;

  /// find and set the corrected abstract item view
  virtual QAbstractItemView* findCorrectedAbstractItemView(QObject* object) const = 0;

protected Q_SLOTS:
  virtual void onClicked(const QModelIndex&);
  virtual void onActivated(const QModelIndex&);
  virtual void onDoubleClicked(const QModelIndex&);
  virtual void onCurrentChanged(const QModelIndex&);
  virtual void onSelectionChanged(const QItemSelection&);

  /// Compute a visual rectangle for the item and signal it
  virtual void onEnteredCheck(const QModelIndex&) = 0;
  virtual void onViewportEnteredCheck();

protected:
  QString getIndexAsString(const QModelIndex& item);
  QString getIndicesAsString(const QModelIndexList& selectedIndices);

  /// Record events from @a abstractItemView when signals are fired.
  /// @sa monitorSignalsInternal()
  void monitorSignals(QAbstractItemView* abstractItemView);

  /// Subclasses (pqTreeViewEventTranslator, possibly others) override this to monitor
  /// additional view specific signals
  virtual void monitorSignalsInternal(QAbstractItemView* abstractItemView) {}

  QPointer<QAbstractItemView> AbstractItemView;
  QPointer<QItemSelectionModel> ItemSelectionModel;
  const QModelIndex* ModelItemCheck;
  bool Editing;
  bool AbstractItemViewMouseTracking;
  bool Checking;
  QPoint LastPos;

private:
  pqAbstractItemViewEventTranslatorBase(
    const pqAbstractItemViewEventTranslatorBase&);              // Not implemented.
  void operator=(const pqAbstractItemViewEventTranslatorBase&); // Not implemented.
};

#endif
