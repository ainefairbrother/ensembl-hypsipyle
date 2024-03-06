"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

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
