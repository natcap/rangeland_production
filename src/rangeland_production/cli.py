# coding=UTF-8
"""Entry point for rangeland production, modeled on InVEST applications."""
from __future__ import absolute_import

import argparse
import os
import importlib
import logging
import sys
import collections
import pprint
import multiprocessing
import json

try:
    from . import __version__
    from . import utils
    from . import datastack
except ImportError:
    # When we're in a pyinstaller build, this isn't a module.
    from rangeland_production import __version__
    from rangeland_production import utils
    from rangeland_production import datastack


DEFAULT_EXIT_CODE = 1
LOGGER = logging.getLogger(__name__)

_PYNAME = 'rangeland_production.forage'
_GUI_CLASS = 'forage.Forage'


def main(user_args=None):
    """CLI entry point for launching rangeland production model.

    This command-line interface supports two methods of launching the rangeland
    production model from the command-line:

        * through its GUI
        * in headless mode, without its GUI.

    Running in headless mode allows us to bypass all GUI functionality,
    so models may be run in this way wthout having GUI packages
    installed.
    """

    parser = argparse.ArgumentParser(description=(
        'Rangeland Production Model.'),
        prog='rangeland_production'
    )
    parser.add_argument('--version', action='version',
                        version=__version__)
    verbosity_group = parser.add_mutually_exclusive_group()

    verbosity_group.add_argument('-v', '--verbose', dest='verbosity', default=0,
                                 action='count', help=(
                                     'Increase verbosity. Affects how much is '
                                     'printed to the console and (if running '
                                     'in headless mode) how much is written '
                                     'to the logfile.'))
    verbosity_group.add_argument('--debug', dest='log_level',
                                 default=logging.CRITICAL,
                                 action='store_const', const=logging.DEBUG,
                                 help='Enable debug logging. Alias for -vvvvv')
    parser.add_argument('-l', '--headless', action='store_true',
                        dest='headless',
                        help=('Attempt to run model without its GUI.'))
    parser.add_argument('-d', '--datastack', default=None, nargs='?',
                        help='Run the specified model with this datastack')
    parser.add_argument('-w', '--workspace', default=None, nargs='?',
                        help='The workspace in which outputs will be saved')

    gui_options_group = parser.add_argument_group(
        'gui options', 'These options are ignored if running in headless mode')
    gui_options_group.add_argument('-q', '--quickrun', action='store_true',
                                   help=('Run the target model without '
                                         'validating and quit with a nonzero '
                                         'exit status if an exception is '
                                         'encountered'))

    cli_options_group = parser.add_argument_group('headless options')
    cli_options_group.add_argument('-y', '--overwrite', action='store_true',
                                   default=False,
                                   help=('Overwrite the workspace without '
                                         'prompting for confirmation'))
    cli_options_group.add_argument('-n', '--no-validate', action='store_true',
                                   dest='validate', default=True,
                                   help=('Do not validate inputs before '
                                         'running the model.'))

    args = parser.parse_args(user_args)

    root_logger = logging.getLogger()
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        fmt='%(asctime)s %(name)-18s %(levelname)-8s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S ')
    handler.setFormatter(formatter)

    # Set the log level based on what the user provides in the available
    # arguments.  Verbosity: the more v's the lower the logging threshold.
    # If --debug is used, the logging threshold is 10.
    # If the user goes lower than logging.DEBUG, default to logging.DEBUG.
    log_level = min(args.log_level, logging.CRITICAL - (args.verbosity*10))
    handler.setLevel(max(log_level, logging.DEBUG))  # don't go lower than DEBUG
    root_logger.addHandler(handler)
    LOGGER.info('Setting handler log level to %s', log_level)

    # FYI: Root logger by default has a level of logging.WARNING.
    # To capture ALL logging produced in this system at runtime, use this:
    # logging.getLogger().setLevel(logging.DEBUG)
    # Also FYI: using logging.DEBUG means that the logger will defer to
    # the setting of the parent logger.
    logging.getLogger('rangeland_production').setLevel(logging.DEBUG)

    # Now that we've set up logging based on args, we can start logging.
    LOGGER.debug(args)

    if args.headless:
        target_mod = _PYNAME
        model_module = importlib.import_module(name=target_mod)
        LOGGER.info('imported target %s from %s',
                    model_module.__name__, model_module)

        paramset = datastack.extract_parameter_set(args.datastack)
        try:
            parsed_datastack = datastack.extract_parameter_set(args.datastack)
        except Exception as error:
            parser.exit(
                1, "Error when parsing JSON datastack:\n    " + str(error))

        if not args.workspace:
            if ('workspace_dir' not in parsed_datastack.args or
                    parsed_datastack.args['workspace_dir'] in ['', None]):
                parser.exit(
                    1, ('Workspace must be defined at the command line '
                        'or in the datastack file'))
        else:
            parsed_datastack.args['workspace_dir'] = args.workspace

        with utils.prepare_workspace(parsed_datastack.args['workspace_dir'],
                                     name=parsed_datastack.model_name,
                                     logging_level=log_level):
            LOGGER.log(datastack.ARGS_LOG_LEVEL,
                       datastack.format_args_dict(paramset.args,
                                                  paramset.model_name))
            if not args.validate:
                LOGGER.info('Skipping validation by user request')
            else:
                model_warnings = []
                try:
                    model_warnings = getattr(
                        model_module, 'validate')(paramset.args)
                except AttributeError:
                    LOGGER.warn(
                        '%s does not have a defined validation function.',
                        paramset.model_name)
                finally:
                    if model_warnings:
                        LOGGER.warn('Warnings found: \n%s',
                                    pprint.pformat(model_warnings))

            if not args.workspace:
                args.workspace = os.getcwd()

            # If the workspace exists and we don't have up-front permission to
            # overwrite the workspace, prompt for permission.
            if (os.path.exists(args.workspace) and
                    len(os.listdir(args.workspace)) > 0 and
                    not args.overwrite):
                overwrite_denied = False
                if not sys.stdout.isatty():
                    overwrite_denied = True
                else:
                    user_response = raw_input(
                        'Workspace exists: %s\n    Overwrite? (y/n) ' % (
                            os.path.abspath(args.workspace)))
                    while user_response not in ('y', 'n'):
                        user_response = raw_input(
                            "Response must be either 'y' or 'n': ")
                    if user_response == 'n':
                        overwrite_denied = True

                if overwrite_denied:
                    # Exit the parser with an error message.
                    parser.exit(DEFAULT_EXIT_CODE,
                                ('Use --workspace to define an '
                                 'alternate workspace.  Aborting.'))
                else:
                    LOGGER.warning(
                        'Overwriting the workspace per user input %s',
                        os.path.abspath(args.workspace))

            if 'workspace_dir' not in paramset.args:
                paramset.args['workspace_dir'] = args.workspace

            # execute the model's execute function with the loaded args
            getattr(model_module, 'execute')(paramset.args)
    else:
        from rangeland_production.ui import inputs
        # import the GUI from the known class
        gui_class = _GUI_CLASS
        module_name, classname = gui_class.split('.')
        module = importlib.import_module(
            name='.ui.%s' % module_name,
            package='rangeland_production')

        # Instantiate the form
        model_form = getattr(module, classname)()

        # load the datastack if one was provided
        try:
            if args.datastack:
                model_form.load_datastack(args.datastack)
        except Exception as error:
            # If we encounter an exception while loading the datastack, log the
            # exception (so it can be seen if we're running with appropriate
            # verbosity) and exit the argparse application with exit code 1 and
            # a helpful error message.
            LOGGER.exception('Could not load datastack')
            parser.exit(DEFAULT_EXIT_CODE,
                        'Could not load datastack: %s\n' % str(error))

        if args.workspace:
            model_form.workspace.set_value(args.workspace)

        # Run the UI's event loop
        model_form.run(quickrun=args.quickrun)
        app_exitcode = inputs.QT_APP.exec_()

        # Handle a graceful exit
        if model_form.form.run_dialog.messageArea.error:
            parser.exit(DEFAULT_EXIT_CODE,
                        'Rangeland Production: run failed\n')

        if app_exitcode != 0:
            parser.exit(app_exitcode,
                        'App terminated with exit code %s\n' % app_exitcode)

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
