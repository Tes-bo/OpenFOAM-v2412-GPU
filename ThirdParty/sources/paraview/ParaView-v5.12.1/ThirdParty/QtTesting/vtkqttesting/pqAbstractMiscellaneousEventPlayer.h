/*=========================================================================

   Program: ParaView
   Module:    pqAbstractMiscellaneousEventPlayer.h

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

#ifndef _pqAbstractMiscellaneousEventPlayer_h
#define _pqAbstractMiscellaneousEventPlayer_h

#include "pqWidgetEventPlayer.h"

class pqTestUtility;

/// Event playback handler for a collection of miscellaneous commands.
/// For these events, the "object" on which the event is triggered is generally
/// immaterial.
///
/// Supported commands are:
/// \li \c pause : pause the application for a fixed time. Time is specified as
///                the arguments in milliseconds to pause the application for.
///                This is a sleep i.e. no events will be processed and the
///                application will appear frozen.
/// \li \c process_events: process pending events (including timers, etc.).
///                Optional time (in milliseconds) may be specified to also pause
///                the test playback. Unlike "pause" however, this will continue
///                to process all events arising in the application e.g.
///                responding to timer events.
/// \li \c dashboard_mode: trigger a call to pqTestUtility::setDashboardMode with the
///                provided argument.

class QTTESTING_EXPORT pqAbstractMiscellaneousEventPlayer : public pqWidgetEventPlayer
{
  Q_OBJECT
  typedef pqWidgetEventPlayer Superclass;

public:
  pqAbstractMiscellaneousEventPlayer(pqTestUtility* util, QObject* p = 0);

  using Superclass::playEvent;
  bool playEvent(
    QObject* Object, const QString& Command, const QString& Arguments, bool& Error) override;

private:
  pqAbstractMiscellaneousEventPlayer(const pqAbstractMiscellaneousEventPlayer&);
  pqAbstractMiscellaneousEventPlayer& operator=(const pqAbstractMiscellaneousEventPlayer&);

  pqTestUtility* TestUtility = nullptr;
};

#endif // !_pqAbstractMiscellaneousEventPlayer_h
