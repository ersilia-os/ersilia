import app1
import app2
import app3
import streamlit as st

PAGES = {"About": app1, "API": app2, "Catalog": app3}
st.sidebar.title("Navigation")
selection = st.sidebar.radio("Go to", list(PAGES.keys()))
page = PAGES[selection]
page.app()
