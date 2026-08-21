---
sidebar_label: core
title: proteolyzer.core
---

The base suite: read search-engine output into a consistent frame.

Submodules
    formats: the input formats recognized, and how their columns are handled
    loader: read a described source into memory
    models: Data, LoadedData and ProcessedData
    processor: dtype narrowing, derived columns, labelling information
    matrix: pivot processed data into a quantitative matrix
    operations: small pure functions
    logging: the package logger and the Logged base class
    io: parquet interchange for frames passed between stages
    pipeline: shared stage plumbing (parameters, progress, provenance)

## frame\_exists

## read\_frame

## write\_frame

## DataLoader

## Logged

## configure\_logging

## MatrixBuilder

## Data

## LoadedData

## ProcessedData

## cv

## NullQueue

## Stage

## DataProcessor

#### \_\_all\_\_

