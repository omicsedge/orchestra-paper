import logging
import os
from typing import List

from . import AncestryModel

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def register_submodel(
    model, sagemaker_session, role, base_image_uri, instance_type, source_root="./"
):
    base_model = SmModel(
        image_uri=base_image_uri,
        model_data=model.base_model_uri,
        sagemaker_session=sagemaker_session,
        role=role,
    )

    smooth_model = PyTorchModel(
        model_data=model.smooth_model_uri,
        sagemaker_session=sagemaker_session,
        role=role,
        framework_version="1.10",
        py_version="py38",
        source_dir=os.path.join(source_root, "smoothing/inference/"),
        entry_point="inference_user.py",
    )

    name = f"ancestry/{model._model.name}/{model._name}".replace(".", "-").replace(
        "/", "--"
    )
    logging.info(f"Model Name: {name}")
    model = PipelineModel(
        name=name,
        role=role,
        sagemaker_session=sagemaker_session,
        models=[base_model, smooth_model],
    )

    model.create(instance_type=instance_type)

    return model


def get_model_registry(model, sub_model_names):
    registry = {
        "artifacts": {
            "population_map_uri": model.population_map_uri,
        },
        "sub-models": {},
    }

    for sm_name in sub_model_names:
        logger.info(f"register {sm_name}")
        submodel = model.get_submodel(sm_name)
        registry["sub-models"][sm_name] = {
            "model": {
                "base_model_uri": submodel.base_model_uri,
                "smooth_model_uri": submodel.smooth_model_uri,
            },
            "artifacts": {
                "parameters": submodel.parameter_file_uri,
                "windows_info": submodel.window_info_file_uri,
            },
        }

    return registry


def register_ancestry_model(model, sub_model_names, *args, **kwargs):
    model_registry = get_model_registry(model, sub_model_names)
    for sm_name, data in model_registry["sub-models"].items():
        logger.info(f"register {sm_name}")
        submodel = model.get_submodel(sm_name)
        sagemaker_model = register_submodel(submodel, *args, **kwargs)
        data["model"]["name"] = sagemaker_model.name

    return model_registry


def register_cfn_models(model: AncestryModel, sub_model_names: List[str]) -> dict:
    submodels = []
    for sm_name in sub_model_names:
        logger.info(f"register {sm_name}")
        submodels.append(model.get_submodel(sm_name))
    return submodels
