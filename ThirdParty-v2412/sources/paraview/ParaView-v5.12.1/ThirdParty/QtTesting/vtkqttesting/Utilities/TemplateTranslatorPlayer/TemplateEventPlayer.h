#ifndef __TemplateEventPlayer_h
#define __TemplateEventPlayer_h

// QtTesting includes
#include <pqWidgetEventPlayer.h>

/// Concrete implementation of pqWidgetEventPlayer that translates
/// high-level events into low-level Qt events.

class TemplateEventPlayer : public pqWidgetEventPlayer
{
  Q_OBJECT
  typedef pqWidgetEventPlayer Superclass;

public:
  TemplateEventPlayer(QObject* parent = 0);

  using Superclass::playEvent;
  bool playEvent(QObject* Object, const QString& Command, const QString& Arguments, bool& Error);

private:
  TemplateEventPlayer(const TemplateEventPlayer&);            // NOT implemented
  TemplateEventPlayer& operator=(const TemplateEventPlayer&); // NOT implemented
};

#endif
