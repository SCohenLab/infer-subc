### Annotated on (25-29)-12-24

from .segmenter_function import SegmenterFunction, FunctionParameter, WidgetType
### SegmenterFunction: Class representing a segmentation function (is this a single workflow step???)
### FunctionParameter: Class representing an input parameter to a segmentation function
### WidgetType: Class establishing whether the widget (input parameter widget) is a slider or dropdown menu

from .workflow_step import WorkflowStep, WorkflowStepCategory
### WorkflowStep: Class representing an individual step in the segmentation workflow
### WorkflowStepCategory: Class establishing the category of the workflow step

from .workflow import Workflow
### Workflow: Class representing an executable segmentation workflow

from .batch_workflow import BatchWorkflow
### BatchWorkflow: Class representing a batch of multiple workflows to process
# from .workflow_definition import WorkflowDefinition (old code)

from .workflow_definition import WorkflowDefinition, PrebuiltWorkflowDefinition
### (this is where the function we want to fix lives)
### WorkflowDefinition: Class that defines a custom workflow loaded from a file
### PrebuiltWorkflowDefinition: Class that defines a prebuilt workflow from using files in the assets folder

from .workflow_engine import WorkflowEngine
### WorkflowEngine: Class used to access and execute aicssegmentation workflows


### There was only one difference, now the PrebuiltWorkflowDefinition class is imported
### All of the imports are from the other py files in the folder (all py files are referenced except for workflow_config)