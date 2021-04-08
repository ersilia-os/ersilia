import streamlit as st
from ..app import AppBase


class StreamlitBlocks(AppBase):
    def __init__(self, model_id, config_json=None):
        AppBase.__init__(self, config_json=config_json)
        self.card = self.get_model_card(model_id)
        self.model_id = model_id

    def load_model(self):
        pass

    def header(self):
        st.title(self.card.title)
        st.markdown("`MODEL ID: %s`" % self.model_id)
        st.markdown(
            self.card.description
            + " "
            + "To know more, please visit [Ersilia Model Hub](%s/%s)."
            % (self.cfg.HUB.WEB, self.model_id)
        )

    @staticmethod
    def molecule_input():
        st.subheader("Input molecule")
        inp = st.text_input(
            "Please enter a molecule (name, SMILES...)",
            value="",
            max_chars=None,
            key=None,
            type="default",
        )
        return inp

    @staticmethod
    def protein_input():
        inp = st.text_input(
            "Please enter a protein (sequence, UniProtAC...)",
            value="",
            max_chars=None,
            key=None,
            type="default",
        )
        return inp
