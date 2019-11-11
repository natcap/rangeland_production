from PyInstaller.utils.hooks import (collect_data_files,
                                     collect_submodules,
                                     copy_metadata)

hiddenimports = collect_submodules('rangeland_production')
datas = copy_metadata(
    "rangeland_production") + collect_data_files('rangeland_production')
