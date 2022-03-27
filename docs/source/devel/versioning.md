# Versioning

**We use semantic versioning.**

- To create a new tag, run the below command:
```
git tag 'v0.0.0' -a -m "Test version"
```
- To push tags to GitHub, run the below command:
```
git push --follow-tags
```

## Publish to PyPi
A GitHub Action is ready to publish to PyPi whenever a release is done. Releases can be created from tags using github web interface.
