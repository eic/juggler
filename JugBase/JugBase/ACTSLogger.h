// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#ifndef RECTRACKER_ACTSLOGGER_H
#define RECTRACKER_ACTSLOGGER_H

#include <Acts/Utilities/Logger.hpp>
#include <GaudiKernel/ServiceHandle.h>
#include <GaudiKernel/IMessageSvc.h>
#include <GaudiKernel/MsgStream.h>

/** Filter
 *
 * \ingroup base
 */
class GaudiFilterPolicy : public Acts::Logging::OutputFilterPolicy {
public:
  GaudiFilterPolicy(IMessageSvc* owner, const Acts::Logging::Level& lvl) : m_owner(owner), m_level(lvl) {}

  bool doPrint(const Acts::Logging::Level& lvl) const {
    
    MSG::Level l = MSG::VERBOSE;
    switch (lvl) {
    case Acts::Logging::VERBOSE:
      l = MSG::VERBOSE;
      break;
    case Acts::Logging::DEBUG:
      l = MSG::DEBUG;
      break;
    case Acts::Logging::INFO:
      l = MSG::INFO;
      break;
    case Acts::Logging::WARNING:
      l = MSG::WARNING;
      break;
    case Acts::Logging::ERROR:
      l = MSG::ERROR;
      break;
    case Acts::Logging::FATAL:
      l = MSG::FATAL;
      break;
     case Acts::Logging::MAX:
      l = MSG::VERBOSE;
      break;
    }

    return l >= m_owner->outputLevel();
  }

  /// Get the level of this filter policy
  /// @return the levele
  Acts::Logging::Level level() const override { return m_level; }

  /// Make a copy of this filter policy with a new level
  /// @param level the new level
  /// @return the new copy
  std::unique_ptr<OutputFilterPolicy> clone(Acts::Logging::Level level) const override {
    return std::make_unique<GaudiFilterPolicy>(m_owner, level);
  }

private:
  IMessageSvc* m_owner;
  Acts::Logging::Level m_level;
};

class GaudiPrintPolicy : public Acts::Logging::OutputPrintPolicy {
public:
  GaudiPrintPolicy(IMessageSvc* owner) : m_owner(owner),m_messenger(owner) {}

  const std::string& name() const override {
    return m_name;
  };

  /// Make a copy of this print policy with a new name
  /// @param name the new name
  /// @return the copy
  std::unique_ptr<Acts::Logging::OutputPrintPolicy> clone(
    const std::string& name) const override {
    (void)name;
    return std::make_unique<GaudiPrintPolicy>(m_owner);
  };

  void setName(std::string name) {
    m_name = name;
  }

  void flush(const Acts::Logging::Level& lvl, const std::string& input) {
    MSG::Level l = MSG::VERBOSE;
    switch (lvl) {
    case Acts::Logging::VERBOSE:
      l = MSG::VERBOSE;
      break;
    case Acts::Logging::DEBUG:
      l = MSG::DEBUG;
      break;
    case Acts::Logging::INFO:
      l = MSG::INFO;
      break;
    case Acts::Logging::WARNING:
      l = MSG::WARNING;
      break;
    case Acts::Logging::ERROR:
      l = MSG::ERROR;
      break;
    case Acts::Logging::FATAL:
      l = MSG::FATAL;
      break;
     case Acts::Logging::MAX:
      l = MSG::VERBOSE;
      break;
    }

    m_messenger << l << m_name << "\t" << input << endmsg;
  }

private:
  IMessageSvc* m_owner;
  MsgStream m_messenger;
  std::string m_name;
};

#endif  // RECTRACKER_ACTSLOGGER_H
