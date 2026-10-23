
from src.dataframe_transformation import filter_dataframe, add_measure_id


class AnalysisInstance:
    """
    Generic pipeline component.

    Designed to standardize the lifecycle of a data-processing step.
    Subclasses should implement the core logic methods.
    """

    def __init__(self, config=None):
        self.config = config or {}
        self.info_path = self.config.get("info_path", None)
        self.adat_path = self.config.get("adat_path", None)
        self.output_data = None

    # -------------------------
    # Data Interface
    # -------------------------

    def load(self):
        """Load input data into the component."""
        self.info_df = pd.read_csv(self.info_path) if self.info_path else None
        self.input_data = None

    def get_output(self):
        """Return processed output."""
        return self.output_data

    # -------------------------
    # Core Logic (override later)
    # -------------------------

    def preprocess(self):
        """Optional preprocessing step."""
        pass

    def process(self):
        """Main processing logic."""
        raise NotImplementedError("Subclasses must implement process()")

    def postprocess(self):
        """Optional cleanup or formatting."""
        pass
    

    # -------------------------
    # Execution
    # -------------------------

    def run(self, data):
        """
        Standard execution pipeline.
        """
        self.load(data)
        self.preprocess()
        self.process()
        self.postprocess()
        return self.get_output()
    
    def concat_adats(self, )

    # -------------------------
    # Debugging
    # -------------------------

    def __repr__(self):
        return f"{self.__class__.__name__}(config={self.config})"