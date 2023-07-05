"""Custom Thoas logger"""

import logging

from pymongo import monitoring


class CommandLogger(monitoring.CommandListener):
    """Logger for MongoDB transactions"""

    def __init__(self, log):
        self.log = log

    def started(self, event):
        self.log.debug(
            "[Request id: %s] Command %s started on server %s",
            event.request_id,
            event.command_name,
            event.connection_id,
        )
        if event.command_name == "find":
            self.log.debug(
                "[Request id: %s] Running query %s",
                event.request_id,
                event.command["filter"],
            )

    def succeeded(self, event):
        self.log.debug(
            "[Request id: %s] Command %s on server %s succeeded in %s microseconds",
            event.request_id,
            event.command_name,
            event.connection_id,
            event.duration_micros,
        )

    def failed(self, event):
        self.log.debug(
            "[Request id: %s] Command %s on server %s failed in %s microseconds",
            event.request_id,
            event.command_name,
            event.connection_id,
            event.duration_micros,
        )
