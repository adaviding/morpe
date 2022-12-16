﻿using System;
using System.Diagnostics.CodeAnalysis;
using Morpe.Validation;

namespace Morpe
{
    /// <summary>
    /// This tells us the current state of the <see cref="CategorizedData"/>, and also encapsulates logic for transitions
    /// through states.
    ///
    /// There are two types of state transitions:
    /// 1.  Spatial conditioning or de-conditioning.
    /// 2.  Polynomial expansion or contraction.
    ///
    /// Note that spatial conditioning always occurs before polynomial expansion.
    /// </summary>
    public class CategorizedDataState
    {
        /// <summary>
        /// If <see cref="IsConditioned"/>, this was the <see cref="SpatialConditioner"/> that was used.
        /// </summary>
        public SpatialConditioner Conditioner;

        /// <summary>
        /// If true, the data has been spatially conditioned, false otherwise.
        /// </summary>
        public bool IsConditioned;

        /// <summary>
        /// If true, the data has gone through a polynomial expansion, false otherwise.
        /// </summary>
        public bool IsExpanded;
        
        /// <summary>
        /// If <see cref="IsExpanded"/>, then this was the polynomial used for expansion.
        /// </summary>
        public Poly Poly;

        public CategorizedDataState Clone()
        {
            CategorizedDataState output = this.MemberwiseClone() as CategorizedDataState;
            output.Poly = output.Poly.Clone();
            output.Conditioner = output.Conditioner.Clone();

            return output;
        }

        /// <summary>
        /// Performs spatial conditioning and updates the state of this instance.  If the instance is not in the correct
        /// state then an exception will be thrown.
        /// </summary>
        /// <param name="data">The data on which spatial conditioning is to be performed.</param>
        public void Condition([NotNull] CategorizedData data)
        {
            Chk.NotNull(this.Conditioner, nameof(this.Conditioner));
            Chk.True(!this.IsConditioned, "The data is already spatially conditioned.");
            Chk.True(!this.IsExpanded, "The polynomial expansion has already occurred.");

            this.Conditioner.Condition(data);
            this.IsConditioned = true;
        }
        
        /// <summary>
        /// Undoes polynomial expansion and updates the state of this instance.  If the instance is not in the correct
        /// state then an exception will be thrown.
        /// </summary>
        /// <param name="data">The data on which spatial conditioning is to be performed.</param>
        public void Contract([NotNull] CategorizedData data)
        {
            Chk.NotNull(this.Poly, nameof(this.Poly));
            Chk.True(this.IsExpanded, "The data is not expanded by polynomial.");
            
            data.Contract();
            this.IsExpanded = false;
        }
        
        /// <summary>
        /// Undoes spatial conditioning and updates the state of this instance.  If the instance is not in the correct
        /// state then an exception will be thrown.
        /// </summary>
        /// <param name="data">The data on which spatial conditioning is to be performed.</param>
        public void Decondition([NotNull] CategorizedData data)
        {
            Chk.NotNull(this.Conditioner, nameof(this.Conditioner));
            Chk.True(this.IsConditioned, "The data is not spatially conditioned.");
            Chk.True(!this.IsExpanded, "The polynomial expansion has already occurred.");

            this.Conditioner.Decondition(data);
            this.IsConditioned = false;
        }

        /// <summary>
        /// Performs polynomial expansion and updates the state of this instance.  If the instance is not in the correct
        /// state then an exception will be thrown.
        /// </summary>
        /// <param name="data">The data on which spatial conditioning is to be performed.</param>
        public void Expand([NotNull] CategorizedData data)
        {
            Chk.NotNull(this.Poly, nameof(this.Poly));
            Chk.True(!this.IsExpanded, "The polynomial expansion has already occurred.");
            
            data.Expand(this.Poly);
            this.IsExpanded = true;
        }
    }
}