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

try:
    from . import utils
except ValueError:
    # When we're in a pyinstaller build, this isn't a module.
    from rangeland_production import utils


DEFAULT_EXIT_CODE = 1
LOGGER = logging.getLogger(__name__)
_UIMETA = collections.namedtuple('UIMeta', 'pyname gui aliases')

_MODEL_UIS = {
    'rangelands': _UIMETA(
        pyname='rangeland_production.forage',
        gui='forage.Forage',
        aliases=('RPM')),
}

# Build up an index mapping aliase to modelname.
# ``modelname`` is the key to the _MODEL_UIS dict, above.
_MODEL_ALIASES = {}
for _modelname, _meta in _MODEL_UIS.items():
    for _alias in _meta.aliases:
        assert _alias not in _MODEL_ALIASES, (
            'Alias %s already defined for model %s') % (
                _alias, _MODEL_ALIASES[_alias])
        _MODEL_ALIASES[_alias] = _modelname


def build_model_list_table():
    """Build a table of model names, aliases and other details.

    This table is a table only in the sense that its contents are aligned
    into columns, but are not separated by a delimited.  This table
    is intended to be printed to stdout.

    Returns:
        A string representation of the formatted table.
    """
    model_names = sorted(_MODEL_UIS.keys())
    max_model_name_length = max(len(name) for name in model_names)
    max_alias_name_length = max(len(', '.join(meta.aliases))
                                for meta in _MODEL_UIS.values())
    template_string = '    {modelname} {aliases}   {usage}'
    strings = ['Available models:']
    for model_name in sorted(_MODEL_UIS.keys()):
        usage_string = '(No GUI available)'
        if _MODEL_UIS[model_name].gui is not None:
            usage_string = ''

        alias_string = ', '.join(_MODEL_UIS[model_name].aliases)
        if alias_string:
            alias_string = '(%s)' % alias_string

        strings.append(template_string.format(
            modelname=model_name.ljust(max_model_name_length),
            aliases=alias_string.ljust(max_alias_name_length),
            usage=usage_string))
    return '\n'.join(strings) + '\n'


class ListModelsAction(argparse.Action):
    """An argparse action to list the available models."""
    def __call__(self, parser, namespace, values, option_string=None):
        """Print the available models and quit the argparse parser.

        See https://docs.python.org/2.7/library/argparse.html#action-classes
        for the full documentation for argparse classes.

        Overridden from argparse.Action.__call__"""
        setattr(namespace, self.dest, self.const)
        parser.exit(message=build_model_list_table())


class SelectModelAction(argparse.Action):
    """Given a possily-ambiguous model string, identify the model to run.

    This is a subclass of ``argparse.Action`` and is executed when the argparse
    interface detects that the user has attempted to select a model by name.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        """Given the user's input, determine which model they're referring to.

        When the user didn't provide a model name, we print the help and exit
        with a nonzero exit code.

        Identifiable model names are:

            * the model name (verbatim) as identified in the keys of _MODEL_UIS
            * a uniquely identifiable prefix for the model name (e.g. "d"
              matches "delineateit", but "fi" matches both "fisheries" and
              "finfish"
            * a known model alias, as registered in _MODEL_UIS

        If no single model can be identified based on these rules, an error
        message is printed and the parser exits with a nonzero exit code.

        See https://docs.python.org/2.7/library/argparse.html#action-classes
        for the full documentation for argparse classes and this __call__
        method.

        Overridden from argparse.Action.__call__"""
        if values in ['', None]:
            parser.print_help()
            parser.exit(1, message=build_model_list_table())
        else:
            known_models = sorted(list(_MODEL_UIS.keys()) + ['launcher'])

            matching_models = [model for model in known_models if
                               model.startswith(values)]

            exact_matches = [model for model in known_models if
                             model == values]

            if len(matching_models) == 1:  # match an identifying substring
                modelname = matching_models[0]
            elif len(exact_matches) == 1:  # match an exact modelname
                modelname = exact_matches[0]
            elif values in _MODEL_ALIASES:  # match an alias
                modelname = _MODEL_ALIASES[values]
            elif len(matching_models) == 0:
                parser.exit(status=1, message=(
                    "Error: '%s' not a known model" % values))
            else:
                parser.exit(
                    status=1,
                    message=(
                        "Model string '{model}' is ambiguous:\n"
                        "    {matching_models}").format(
                            model=values,
                            matching_models=' '.join(matching_models)))
        setattr(namespace, self.dest, modelname)


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
    list_group = parser.add_mutually_exclusive_group()
    verbosity_group = parser.add_mutually_exclusive_group()
    import rangeland_production

    parser.add_argument('--version', action='version',
                        version='0.1.0')  # TODO rangeland_production.__version__)
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
    list_group.add_argument('--list', action=ListModelsAction,
                            nargs=0, const=True,
                            help='List available models')
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

    list_group.add_argument('model', action=SelectModelAction, nargs='?',
                            help=('The model/tool to run. Use --list to show '
                                  'available models/tools. Identifiable model '
                                  'prefixes may also be used. Alternatively,'
                                  'specify "launcher" to reveal a model '
                                  'launcher window.'))

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

    try:
        # Importing model UI files here will usually import qtpy before we can
        # set the sip API in natcap.invest.ui.inputs.
        # Set it here, before we can do the actual importing.
        import sip
        # 2 indicates SIP/Qt API version 2
        sip.setapi('QString', 2)

        from rangeland_production.ui import inputs
    except ImportError as error:
        # Can't import UI, exit with nonzero exit code
        LOGGER.exception('Unable to import the UI')
        parser.error(('Unable to import the UI (failed with "%s")\n'
                      'Is the UI installed?\n'
                      '    pip install rangeland_production[ui]') % error)

    if args.model == 'launcher':
        from rangeland_production.ui import launcher
        launcher.main()

    elif args.headless:
        from rangeland_production import datastack
        target_mod = _MODEL_UIS[args.model].pyname
        model_module = importlib.import_module(name=target_mod)
        LOGGER.info('imported target %s from %s',
                    model_module.__name__, model_module)

        paramset = datastack.extract_parameter_set(args.datastack)

        # prefer CLI option for workspace dir, but use paramset workspace if
        # the CLI options do not define a workspace.
        if args.workspace:
            workspace = os.path.abspath(args.workspace)
            paramset.args['workspace_dir'] = workspace
        else:
            if 'workspace_dir' in paramset.args:
                workspace = paramset.args['workspace_dir']
            else:
                parser.exit(DEFAULT_EXIT_CODE, (
                    'Workspace not defined. \n'
                    'Use --workspace to specify or add a '
                    '"workspace_dir" parameter to your datastack.'))

        with utils.prepare_workspace(workspace,
                                     name=paramset.model_name,
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
        # import the GUI from the known class
        gui_class = _MODEL_UIS[args.model].gui
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
                        'Model %s: run failed\n' % args.model)

        if app_exitcode != 0:
            parser.exit(app_exitcode,
                        'App terminated with exit code %s\n' % app_exitcode)

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()